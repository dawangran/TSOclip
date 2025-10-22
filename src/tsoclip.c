// TSOclip: fast TSO detection & trimming in long reads
// - Tail-window search (default 150bp)
// - PARTIAL = TSO PREFIX anywhere within tail window
// - FULL with ±shift small indels
// - OMP parallel; gzip I/O; single-write per read; optional plain text out
//
// Build: gcc -O3 -march=native -pipe -fopenmp -DNDEBUG -o tsoclip src/tsoclip.c -lz
// License: MIT

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>
#include <errno.h>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/*************** minimal kseq for gz FASTQ ***************/
typedef struct { unsigned char *buf; int begin, end, is_eof; gzFile f; } kstream_t;
typedef struct { size_t l, m; char *s; } kstring_t;
typedef struct { kstring_t name, comment, seq, qual; int last_char; kstream_t *f; } kseq_t;
static inline kstream_t *ks_init(gzFile f){ kstream_t *ks=(kstream_t*)calloc(1,sizeof(kstream_t)); ks->f=f; ks->buf=(unsigned char*)malloc(1<<18); ks->begin=ks->end=0; ks->is_eof=0; return ks; }
static inline void ks_destroy(kstream_t *ks){ if(!ks) return; free(ks->buf); free(ks); }
static inline int ks_getc(kstream_t *ks){
    if (ks->begin >= ks->end){
        if (ks->is_eof) return -1;
        ks->begin=0; ks->end=gzread(ks->f, ks->buf, 1<<18);
        if (ks->end < (1<<18)) ks->is_eof = 1;
        if (ks->end==0) return -1;
    }
    return (int)ks->buf[ks->begin++];
}
static inline kseq_t *kseq_init(gzFile fd){ kseq_t *s=(kseq_t*)calloc(1,sizeof(kseq_t)); s->f=ks_init(fd); s->last_char=0; return s; }
static inline void kseq_destroy(kseq_t *ks){ if(!ks) return; free(ks->name.s); free(ks->comment.s); free(ks->seq.s); free(ks->qual.s); ks_destroy(ks->f); free(ks); }
static inline int kseq_read(kseq_t *seq){
    int c; kstring_t *s=&seq->seq, *q=&seq->qual, *n=&seq->name, *cm=&seq->comment;
    if (seq->last_char == 0){ while ((c=ks_getc(seq->f))!=-1 && c!='@'); if (c==-1) return -1; seq->last_char=c; }
    n->l=0; if ((c=ks_getc(seq->f))==-1) return -1;
    while(c!=-1 && c!='\n' && c!='\r' && c!=' '){ if (n->l+1>=n->m){ n->m=n->m? n->m<<1:256; n->s=(char*)realloc(n->s,n->m);} n->s[n->l++]=c; c=ks_getc(seq->f); }
    n->s=(char*)realloc(n->s,n->l+1); n->s[n->l]=0; while(c!=-1 && c!='\n') c=ks_getc(seq->f);
    s->l=0; while ((c=ks_getc(seq->f))!=-1 && c!='+' && c!='@'){ if (c=='\n'||c=='\r') continue; if (s->l+1>=s->m){ s->m=s->m? s->m<<1:1024; s->s=(char*)realloc(s->s,s->m);} s->s[s->l++]=c; }
    if (c=='+'){ while ((c=ks_getc(seq->f))!=-1 && c!='\n'); q->l=0; while (q->l < s->l && (c=ks_getc(seq->f))!=-1){
            if (c=='\n'||c=='\r') continue; if (q->l+1>=q->m){ q->m=q->m? q->m<<1:1024; q->s=(char*)realloc(q->s,q->m);} q->s[q->l++]=c; }
        if (q->l != s->l){ q->s=(char*)realloc(q->s, s->l+1); while (q->l < s->l) q->s[q->l++]='I'; }
        q->s[q->l]=0; seq->last_char=0; return (int)s->l;
    } else { seq->last_char=c; return (int)s->l; }
}

/*************** utils ***************/
static inline void to_upper(char *s, int l){ for(int i=0;i<l;i++){ unsigned char c=s[i]; if (c>='a' && c<='z') s[i]=c-32; } }

/*************** TSO-only scanner (tail -> head) ***************/
typedef struct { int start, end, mm, overlap; double mmr; int partial; } TSOHit;
typedef struct { TSOHit *a; int n, m; } TSOHits;
static void th_init(TSOHits *v){ v->a=NULL; v->n=v->m=0; }
static void th_free(TSOHits *v){ free(v->a); v->a=NULL; v->n=v->m=0; }
static void th_push(TSOHits *v, TSOHit h){ if(v->n==v->m){ v->m=v->m? v->m<<1:64; v->a=(TSOHit*)realloc(v->a,v->m*sizeof(TSOHit)); } v->a[v->n++]=h; }

typedef struct {
    int tail_window, tso_max_mm, tso_max_shift, tso_min_overlap, tso_max_hits, min_spacing;
    double tso_max_mmr;
    int n_as_match;
} TSOParams;

static inline int hamming_mm(const char *a, const char *b, int L, int n_as_match){
    int mm=0;
    if (!n_as_match){ for(int i=0;i<L;i++) if(a[i]!=b[i]) mm++; }
    else{ for(int i=0;i<L;i++){ char qa=a[i], qb=b[i]; if (qa=='N') continue; if (qa!=qb) mm++; } }
    return mm;
}

static int cmp_desc(const void *A, const void *B){
    const TSOHit *a=(const TSOHit*)A, *b=(const TSOHit*)B;
    if (a->start!=b->start) return (a->start>b->start)?-1:1;
    if (a->end!=b->end)     return (a->end>b->end)?-1:1;
    if (a->partial!=b->partial) return (a->partial<b->partial)?-1:1;
    if (a->mmr!=b->mmr)     return (a->mmr<b->mmr)?-1:1;
    return 0;
}

/*** PARTIAL: TSO PREFIX anywhere in tail window ***/
static void partial_prefix_hits(const char *seq, int L, const char *tso, int tL,
                                const TSOParams *P, TSOHits *hv){
    int win_start = (L - P->tail_window > 0) ? L - P->tail_window : 0;
    int win_len   = L - win_start;
    for (int k = P->tso_min_overlap; k <= tL - 1; k++){
        const char *tso_pre = tso;
        for (int s = 0; s + k <= win_len; s++){
            int mm = hamming_mm(seq + win_start + s, tso_pre, k, P->n_as_match);
            double mmr = (k > 0) ? (double)mm / k : 1.0;
            if (mmr <= P->tso_max_mmr){
                TSOHit h = { win_start + s, win_start + s + k, mm, k, mmr, 1 };
                th_push(hv, h);
            }
        }
    }
}

/*** FULL: allow ±shift small indels (scan in tail window) ***/
static void full_hits(const char *seq, int L, const char *tso, int tL, const TSOParams *P, TSOHits *hv){
    int win_start = (L - P->tail_window > 0) ? L - P->tail_window : 0, win_len = L - win_start;
    for(int i=0; i<=win_len - tL; i++){
        int best_mm = 1e9;
        for(int sh=-P->tso_max_shift; sh<=P->tso_max_shift; sh++){
            int a = i + (sh<0 ? sh : 0), b = a + tL; if (a<0 || b>win_len) continue;
            int mm = hamming_mm(seq + win_start + a, tso, tL, P->n_as_match);
            if (mm < best_mm) best_mm = mm;
        }
        if (best_mm <= P->tso_max_mm){ int start = win_start + i, end = start + tL; TSOHit h={start,end,best_mm,tL,(double)best_mm/tL,0}; th_push(hv,h); }
    }
}

static void dedup_prune(TSOHits *hv, int min_spacing, int max_keep){
    if (hv->n==0) return; qsort(hv->a, hv->n, sizeof(TSOHit), cmp_desc);
    int w=0; for(int i=0;i<hv->n;i++){
        if(w==0){ hv->a[w++]=hv->a[i]; continue; }
        TSOHit *pr=&hv->a[w-1], *cu=&hv->a[i];
        if (pr->start - cu->start < min_spacing){
            int keep_prev = (pr->mmr < cu->mmr) || (pr->mmr==cu->mmr && pr->mm <= cu->mm);
            if (!keep_prev) hv->a[w-1] = *cu;
        }else hv->a[w++] = *cu;
    }
    hv->n = w; if (max_keep>0 && hv->n > max_keep) hv->n = max_keep;
}

/*** NEW: choose earliest cut-start to remove ALL concatenated TSOs in window ***/
static inline int earliest_cut_start(const TSOHits *hv){
    // hv->n > 0 guaranteed by caller
    int cut_start = hv->a[0].start;
    int have_full = 0;
    for (int k = 0; k < hv->n; ++k){
        const TSOHit *h = &hv->a[k];
        if (h->partial == 0){ // FULL takes priority
            if (!have_full || h->start < cut_start) cut_start = h->start;
            have_full = 1;
        } else if (!have_full && h->start < cut_start){
            // Only consider PARTIAL for earliest cut when no FULL exists
            cut_start = h->start;
        }
    }
    return cut_start;
}

static void tso_scan(const char *seq, int L, const char *tso, int tL, const TSOParams *P, TSOHits *hv){
    th_init(hv); partial_prefix_hits(seq,L,tso,tL,P,hv); full_hits(seq,L,tso,tL,P,hv); dedup_prune(hv,P->min_spacing,P->tso_max_hits);
}

/*************** CLI & options ***************/
typedef struct {
    const char *fastq, *tso, *out_tsv, *out_trim_fastq;
    int tail_window, tso_max_mm, tso_max_shift, tso_min_overlap, tso_max_hits, min_spacing;
    double tso_max_mmr;
    int threads, batch_size;
    int concat_prefer_full;
    int emit_only_hit, min_keep_len, n_as_match;
    int no_json;
    int gzbuf_kb;
    int gzip_level;
    int plain_out;
} Opt;

static void usage(){
    fprintf(stderr,
"TSOclip v0.1.0-mod\n"
"Usage:\n"
"  tsoclip --fastq in.fastq[.gz|-/stdin] --tso SEQ --out-tsv hits.tsv --out-trim-fastq out.fastq[.gz|-/stdout]\n"
"Scan:\n"
"  --tail-window INT      [150]\n"
"  --tso-max-mm INT       [5]\n"
"  --tso-max-shift INT    [4]\n"
"  --tso-min-overlap INT  [12]\n"
"  --tso-max-mmr FLOAT    [0.20]\n"
"  --tso-max-hits INT     [100]\n"
"  --min-spacing INT      [6]\n"
"Behavior:\n"
"  --n-as-match                 # read 'N' as wildcard\n"
"  --concat-prefer-full         # prefer nearest FULL when concatenated & tail-most is partial (reporting only)\n"
"  --emit-only-hit 0|1          # write only reads with accepted hit [0]\n"
"  --min-keep-len INT           # drop reads shorter than this after trimming [1]\n"
"Output speed:\n"
"  --no-json                    # skip all_hits_json (faster)\n"
"  --plain-out                  # write plain text FASTQ (then pipe to pigz yourself)\n"
"  --gzbuf-kb INT               # zlib buffer KB [1024]\n"
"  --gzip-level INT             # gzip level [1]\n"
"Parallel:\n"
"  --threads INT [auto], --batch-size INT [40000]\n");
}
static int parse_int(const char *s){ return (int)strtol(s,NULL,10); }
static double parse_d(const char *s){ return strtod(s,NULL); }

static int parse_args(int argc, char **argv, Opt *o){
    memset(o,0,sizeof(*o));
    o->tail_window=150; o->tso_max_mm=5; o->tso_max_shift=4; o->tso_min_overlap=12;
    o->tso_max_mmr=0.20; o->tso_max_hits=100; o->min_spacing=6;
    o->threads=0; o->batch_size=40000;
    o->emit_only_hit=0; o->min_keep_len=1; o->n_as_match=1; o->concat_prefer_full=0;
    o->no_json=0; o->gzbuf_kb=1024; o->gzip_level=1; o->plain_out=0;
    for(int i=1;i<argc;i++){
        if (!strcmp(argv[i],"--fastq") && i+1<argc) o->fastq=argv[++i];
        else if(!strcmp(argv[i],"--tso") && i+1<argc) o->tso=argv[++i];
        else if(!strcmp(argv[i],"--out-tsv") && i+1<argc) o->out_tsv=argv[++i];
        else if(!strcmp(argv[i],"--out-trim-fastq") && i+1<argc) o->out_trim_fastq=argv[++i];
        else if(!strcmp(argv[i],"--tail-window") && i+1<argc) o->tail_window=parse_int(argv[++i]);
        else if(!strcmp(argv[i],"--tso-max-mm") && i+1<argc) o->tso_max_mm=parse_int(argv[++i]);
        else if(!strcmp(argv[i],"--tso-max-shift") && i+1<argc) o->tso_max_shift=parse_int(argv[++i]);
        else if(!strcmp(argv[i],"--tso-min-overlap") && i+1<argc) o->tso_min_overlap=parse_int(argv[++i]);
        else if(!strcmp(argv[i],"--tso-max-mmr") && i+1<argc) o->tso_max_mmr=parse_d(argv[++i]);
        else if(!strcmp(argv[i],"--tso-max-hits") && i+1<argc) o->tso_max_hits=parse_int(argv[++i]);
        else if(!strcmp(argv[i],"--min-spacing") && i+1<argc) o->min_spacing=parse_int(argv[++i]);
        else if(!strcmp(argv[i],"--threads") && i+1<argc) o->threads=parse_int(argv[++i]);
        else if(!strcmp(argv[i],"--batch-size") && i+1<argc) o->batch_size=parse_int(argv[++i]);
        else if(!strcmp(argv[i],"--concat-prefer-full")) o->concat_prefer_full=1;
        else if(!strcmp(argv[i],"--n-as-match")) o->n_as_match=1;
        else if(!strcmp(argv[i],"--emit-only-hit") && i+1<argc) o->emit_only_hit=parse_int(argv[++i]);
        else if(!strcmp(argv[i],"--min-keep-len") && i+1<argc) o->min_keep_len=parse_int(argv[++i]);
        else if(!strcmp(argv[i],"--no-json")) o->no_json=1;
        else if(!strcmp(argv[i],"--gzbuf-kb") && i+1<argc) o->gzbuf_kb=parse_int(argv[++i]);
        else if(!strcmp(argv[i],"--gzip-level") && i+1<argc) o->gzip_level=parse_int(argv[++i]);
        else if(!strcmp(argv[i],"--plain-out")) o->plain_out=1;
        else { usage(); return -1; }
    }
    if (!o->fastq || !o->tso || !o->out_tsv || !o->out_trim_fastq){ usage(); return -1; }
    if (o->gzip_level<1 || o->gzip_level>9) o->gzip_level=1;
    return 0;
}

/*************** data structs + reusable buffers ***************/
typedef struct { char *id, *seq, *qual; int len; } ReadRec;
typedef struct { int hit_count; TSOHit pick; char *all_json; int use_hit; int out_len; } ReadRes;

typedef struct { char *buf; int cap; } DynBuf;
static inline void db_init(DynBuf *b){ b->buf=NULL; b->cap=0; }
static inline void db_ens(DynBuf *b, int need){
    if (b->cap < need){ int nc = b->cap? b->cap:1024; while (nc < need) nc <<= 1; b->buf=(char*)realloc(b->buf, nc); b->cap=nc; }
}
static inline void db_free(DynBuf *b){ free(b->buf); b->buf=NULL; b->cap=0; }

static void free_batch(ReadRec *R, ReadRes *S, int n, int reuse){
    if (R){
        for(int i=0;i<n;i++){
            free(R[i].id);
            free(R[i].seq);
            free(R[i].qual);
        }
        if(!reuse) free(R);
    }
    if (S){
        for(int i=0;i<n;i++) free(S[i].all_json);
        if(!reuse) free(S);
    }
}

/*************** main ***************/
int main(int argc, char **argv){
    Opt opt; if (parse_args(argc, argv, &opt)!=0) return 1;

    int tL = (int)strlen(opt.tso);
    char *TSO = strdup(opt.tso); for(char *p=TSO; *p; ++p) *p=toupper((unsigned char)*p);
    TSOParams P = { opt.tail_window, opt.tso_max_mm, opt.tso_max_shift,
                    opt.tso_min_overlap, opt.tso_max_hits, opt.min_spacing,
                    opt.tso_max_mmr, opt.n_as_match };

#ifdef _OPENMP
    if (opt.threads>0) omp_set_num_threads(opt.threads);
#endif

    gzFile fp = NULL;
    if (!strcmp(opt.fastq,"-") || !strcmp(opt.fastq,"/dev/stdin")) fp = gzdopen(fileno(stdin), "rb");
    else fp = gzopen(opt.fastq, "rb");
    if(!fp){ fprintf(stderr,"[ERR] open %s: %s\n", opt.fastq, strerror(errno)); return 2; }
    if (opt.gzbuf_kb>0) gzbuffer(fp, opt.gzbuf_kb<<10);

    kseq_t *ks = kseq_init(fp);

    FILE *tsv = fopen(opt.out_tsv, "w"); if(!tsv){ fprintf(stderr,"[ERR] write %s\n", opt.out_tsv); return 3; }

    gzFile fq_gz = NULL; FILE *fq_plain = NULL; int use_gz = !opt.plain_out;
    if (use_gz){
        if (!strcmp(opt.out_trim_fastq,"-") || !strcmp(opt.out_trim_fastq,"/dev/stdout")) fq_gz = gzdopen(fileno(stdout), "wb");
        else fq_gz = gzopen(opt.out_trim_fastq, "wb");
        if(!fq_gz){ fprintf(stderr,"[ERR] write %s\n", opt.out_trim_fastq); return 4; }
        if (opt.gzbuf_kb>0) gzbuffer(fq_gz, opt.gzbuf_kb<<10);
        gzsetparams(fq_gz, opt.gzip_level, Z_DEFAULT_STRATEGY);
    } else {
        if (!strcmp(opt.out_trim_fastq,"-") || !strcmp(opt.out_trim_fastq,"/dev/stdout")) fq_plain = stdout;
        else fq_plain = fopen(opt.out_trim_fastq, "wb");
        if(!fq_plain){ fprintf(stderr,"[ERR] write %s\n", opt.out_trim_fastq); return 4; }
    }

    fprintf(tsv, "read_id\tread_len\ttso_hit_count\tpick_start\tpick_end\tpick_mm\tpick_overlap\tpick_mmr\tpick_partial\tall_hits_json\n");

    const int B = opt.batch_size;
    ReadRec *R = (ReadRec*)calloc(B, sizeof(ReadRec));
    ReadRes *S = (ReadRes*)calloc(B, sizeof(ReadRes));
    DynBuf wbuf; db_init(&wbuf);

    long long total=0, written=0;

    for(;;){
        int n=0;
        for(; n<B; ){
            int r = kseq_read(ks);
            if (r < 0) break;
            R[n].id  = strdup(ks->name.s);
            R[n].seq = strndup(ks->seq.s, ks->seq.l);
            R[n].qual= strndup(ks->qual.s, ks->qual.l);
            R[n].len = (int)ks->seq.l;
            to_upper(R[n].seq, R[n].len);
            n++;
        }
        if (n==0) break;

        #pragma omp parallel for schedule(static)
        for(int i=0;i<n;i++){
            TSOHits hv; tso_scan(R[i].seq, R[i].len, TSO, tL, &P, &hv);
            S[i].hit_count = hv.n; S[i].out_len = R[i].len; S[i].all_json = NULL; S[i].use_hit = 0;

            if (hv.n>0){
                TSOHit chosen = hv.a[0];
                if (opt.concat_prefer_full && hv.n>=2 && chosen.partial){
                    for (int k=0;k<hv.n;k++){ if (hv.a[k].partial==0){ chosen = hv.a[k]; break; } }
                }
                S[i].use_hit=1; S[i].pick=chosen;

                /* MODIFIED: cut to the earliest start among hits in window
                   (prefer FULL; fallback to earliest PARTIAL if no FULL)
                   -> remove ALL concatenated TSOs at once */
                S[i].out_len = earliest_cut_start(&hv);

                if (!opt.no_json){
                    int est = hv.n*40 + 16; S[i].all_json = (char*)malloc(est);
                    char *w = S[i].all_json; int rem=est;
                    w += snprintf(w, rem, "["); rem = est - (w - S[i].all_json);
                    for(int k=0;k<hv.n;k++){
                        w += snprintf(w, rem, "%s{\"s\":%d,\"e\":%d,\"mm\":%d,\"ov\":%d,\"p\":%d}",
                                      (k?",":""), hv.a[k].start, hv.a[k].end, hv.a[k].mm, hv.a[k].overlap, hv.a[k].partial);
                        rem = est - (w - S[i].all_json);
                        if (rem < 64){ int used=(int)(w - S[i].all_json); est<<=1; S[i].all_json=(char*)realloc(S[i].all_json, est); w=S[i].all_json+used; rem=est-used; }
                    }
                    snprintf(w, rem, "]");
                } else {
                    S[i].all_json = strdup("[]");
                }
            } else {
                S[i].pick.start=S[i].pick.end=S[i].pick.mm=S[i].pick.overlap=0; S[i].pick.mmr=0.0; S[i].pick.partial=0;
                S[i].all_json = opt.no_json ? strdup("[]") : strdup("[]");
            }
            th_free(&hv);
        }

        for(int i=0;i<n;i++){
            if (S[i].hit_count>0 && S[i].use_hit){
                fprintf(tsv, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.3f\t%s\t%s\n",
                    R[i].id, R[i].len, S[i].hit_count,
                    S[i].pick.start, S[i].pick.end, S[i].pick.mm, S[i].pick.overlap,
                    S[i].pick.mmr, S[i].pick.partial? "YES":"NO", S[i].all_json);
            } else if (S[i].hit_count>0 && !S[i].use_hit){
                fprintf(tsv, "%s\t%d\t%d\t-\t-\t-\t-\t-\t-\t%s\n",
                    R[i].id, R[i].len, S[i].hit_count, S[i].all_json);
            } else {
                fprintf(tsv, "%s\t%d\t0\t-\t-\t-\t-\t-\t-\t[]\n", R[i].id, R[i].len);
            }

            int write_ok = 1;
            if (opt.emit_only_hit && !S[i].use_hit) write_ok = 0;
            if (S[i].out_len < opt.min_keep_len) write_ok = 0;
            if (write_ok){
                int L = S[i].out_len, idlen = (int)strlen(R[i].id);
                int need = 1 + idlen + 1 + L + 1 + 1 + 1 + L + 1;
                db_ens(&wbuf, need);
                char *p = wbuf.buf;
                *p++='@'; memcpy(p, R[i].id, idlen); p+=idlen; *p++='\n';
                memcpy(p, R[i].seq, L); p+=L; *p++='\n';
                *p++='+'; *p++='\n';
                memcpy(p, R[i].qual, L); p+=L; *p++='\n';
                if (use_gz) gzwrite(fq_gz, wbuf.buf, (unsigned)(p - wbuf.buf));
                else fwrite(wbuf.buf, 1, (size_t)(p - wbuf.buf), fq_plain);
                written++;
            }
        }

        total += n;
        free_batch(R,S,n,1);
        memset(R, 0, sizeof(ReadRec)*B);
        memset(S, 0, sizeof(ReadRes)*B);
    }

    free_batch(R,S,0,0); db_free(&wbuf);
    kseq_destroy(ks); gzclose(fp);
    fclose(tsv);
    if (use_gz) gzclose(fq_gz); else if (fq_plain && fq_plain!=stdout) fclose(fq_plain);
    free(TSO);

    fprintf(stderr, "[STAT] reads=%lld, trimmed_written=%lld\n", total, written);
    return 0;
}
