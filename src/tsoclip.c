// TSOclip: fast TSO detection & trimming in long reads
// - Tail-window search (default 150bp)
// - PARTIAL = TSO PREFIX anywhere within tail window
// - FULL/PARTIAL hits scored with bounded edit distance
// - Optional OMP parallel; gzip I/O; single-write per read; optional plain text out
//
// Build: gcc -O3 -march=native -pipe -fopenmp -DNDEBUG -o tsoclip src/tsoclip.c -lz
// License: MIT

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <zlib.h>
#include <errno.h>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef TSOCLIP_VERSION
#define TSOCLIP_VERSION "0.2.0-dev"
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
    int c; kstring_t *s=&seq->seq, *q=&seq->qual, *n=&seq->name;
    if (seq->last_char == 0){ while ((c=ks_getc(seq->f))!=-1 && c!='@'){} if (c==-1) return -1; seq->last_char=c; }
    n->l=0; if ((c=ks_getc(seq->f))==-1) return -1;
    while(c!=-1 && c!='\n' && c!='\r' && c!=' '){ if (n->l+1>=n->m){ n->m=n->m? n->m<<1:256; n->s=(char*)realloc(n->s,n->m);} n->s[n->l++]=c; c=ks_getc(seq->f); }
    n->s=(char*)realloc(n->s,n->l+1); n->s[n->l]=0; while(c!=-1 && c!='\n') c=ks_getc(seq->f);
    s->l=0; while ((c=ks_getc(seq->f))!=-1 && c!='+' && c!='@'){ if (c=='\n'||c=='\r') continue; if (s->l+1>=s->m){ s->m=s->m? s->m<<1:1024; s->s=(char*)realloc(s->s,s->m);} s->s[s->l++]=c; }
    if (c=='+'){ while ((c=ks_getc(seq->f))!=-1 && c!='\n'){} q->l=0; while (q->l < s->l && (c=ks_getc(seq->f))!=-1){
            if (c=='\n'||c=='\r') continue; if (q->l+1>=q->m){ q->m=q->m? q->m<<1:1024; q->s=(char*)realloc(q->s,q->m);} q->s[q->l++]=c; }
        if (q->l != s->l){ q->s=(char*)realloc(q->s, s->l+1); while (q->l < s->l) q->s[q->l++]='I'; }
        q->s[q->l]=0; seq->last_char=0; return (int)s->l;
    } else { seq->last_char=c; return (int)s->l; }
}

/*************** utils ***************/
static inline void to_upper(char *s, int l){ for(int i=0;i<l;i++){ unsigned char c=s[i]; if (c>='a' && c<='z') s[i]=c-32; } }

/*************** TSO-only scanner (tail -> head) ***************/
typedef struct { int edits, mm, ins, del; } AlignScore;
typedef struct {
    int start, end;
    int edits, mm, ins, del;
    int overlap;
    double err_rate;
    int partial;
} TSOHit;
typedef struct { TSOHit *a; int n, m; } TSOHits;
static void th_init(TSOHits *v){ v->a=NULL; v->n=v->m=0; }
static void th_free(TSOHits *v){ free(v->a); v->a=NULL; v->n=v->m=0; }
static void th_push(TSOHits *v, TSOHit h){ if(v->n==v->m){ v->m=v->m? v->m<<1:64; v->a=(TSOHit*)realloc(v->a,v->m*sizeof(TSOHit)); } v->a[v->n++]=h; }

typedef struct {
    int tail_window, tso_max_mm, tso_max_shift, tso_min_overlap, tso_max_hits, min_spacing;
    int max_tail_gap, max_chain_gap;
    double tso_max_mmr;
    int n_as_match;
} TSOParams;

static inline int imin(int a, int b){ return a < b ? a : b; }
static inline int imax(int a, int b){ return a > b ? a : b; }

static inline AlignScore score_inf(void){
    AlignScore s = { INT_MAX/4, INT_MAX/4, INT_MAX/4, INT_MAX/4 };
    return s;
}

static inline int score_is_inf(AlignScore s){ return s.edits >= INT_MAX/8; }

static inline int score_better(AlignScore a, AlignScore b){
    if (a.edits != b.edits) return a.edits < b.edits;
    if ((a.ins + a.del) != (b.ins + b.del)) return (a.ins + a.del) < (b.ins + b.del);
    if (a.mm != b.mm) return a.mm < b.mm;
    if (a.ins != b.ins) return a.ins < b.ins;
    return a.del < b.del;
}

static inline int base_equal(char a, char b, int n_as_match){
    if (n_as_match && (a=='N' || b=='N')) return 1;
    return a == b;
}

static inline AlignScore score_sub(AlignScore s, int mismatch){
    if (score_is_inf(s)) return s;
    if (mismatch){ s.edits++; s.mm++; }
    return s;
}

static inline AlignScore score_ins(AlignScore s){
    if (score_is_inf(s)) return s;
    s.edits++; s.ins++;
    return s;
}

static inline AlignScore score_del(AlignScore s){
    if (score_is_inf(s)) return s;
    s.edits++; s.del++;
    return s;
}

static AlignScore banded_edit_score(const char *query, int qL, const char *ref, int rL,
                                    int band, int n_as_match){
    enum { ALIGN_STACK_MAX = 512 };
    AlignScore prev_stack[ALIGN_STACK_MAX + 1], cur_stack[ALIGN_STACK_MAX + 1];
    AlignScore *prev = prev_stack, *cur = cur_stack;

    if (rL > ALIGN_STACK_MAX){
        prev = (AlignScore*)malloc((size_t)(rL + 1) * sizeof(AlignScore));
        cur = (AlignScore*)malloc((size_t)(rL + 1) * sizeof(AlignScore));
        if (!prev || !cur){
            free(prev); free(cur);
            return score_inf();
        }
    }

    if (band < 0) band = imax(qL, rL);
    for (int j=0; j<=rL; ++j) prev[j] = score_inf();
    prev[0] = (AlignScore){0,0,0,0};
    for (int j=1; j<=imin(rL, band); ++j) prev[j] = score_del(prev[j-1]);

    for (int i=1; i<=qL; ++i){
        for (int j=0; j<=rL; ++j) cur[j] = score_inf();
        int jlo = imax(0, i - band);
        int jhi = imin(rL, i + band);

        if (jlo == 0) cur[0] = score_ins(prev[0]);
        for (int j=imax(1, jlo); j<=jhi; ++j){
            AlignScore best = score_inf();
            AlignScore diag = score_sub(prev[j-1], !base_equal(query[i-1], ref[j-1], n_as_match));
            AlignScore up = score_ins(prev[j]);
            AlignScore left = score_del(cur[j-1]);

            if (score_better(diag, best)) best = diag;
            if (score_better(up, best)) best = up;
            if (score_better(left, best)) best = left;
            cur[j] = best;
        }

        AlignScore *tmp = prev; prev = cur; cur = tmp;
    }

    AlignScore ans = prev[rL];
    if (rL > ALIGN_STACK_MAX){ free(prev); free(cur); }
    return ans;
}

static int cmp_desc(const void *A, const void *B){
    const TSOHit *a=(const TSOHit*)A, *b=(const TSOHit*)B;
    if (a->start!=b->start) return (a->start>b->start)?-1:1;
    if (a->end!=b->end)     return (a->end>b->end)?-1:1;
    if (a->partial!=b->partial) return (a->partial<b->partial)?-1:1;
    if (a->err_rate!=b->err_rate) return (a->err_rate<b->err_rate)?-1:1;
    if (a->edits!=b->edits) return (a->edits<b->edits)?-1:1;
    return 0;
}

/*** PARTIAL: TSO PREFIX anywhere in tail window ***/
static void partial_prefix_hits(const char *seq, int L, const char *tso, int tL,
                                const TSOParams *P, TSOHits *hv){
    int win_start = (L - P->tail_window > 0) ? L - P->tail_window : 0;
    int win_len   = L - win_start;
    for (int k = P->tso_min_overlap; k <= tL - 1; k++){
        const char *tso_pre = tso;
        int max_edits = (int)(P->tso_max_mmr * k + 1e-9);
        int q_min = imax(1, k - P->tso_max_shift);
        int q_max = k + P->tso_max_shift;
        for (int s = 0; s < win_len; s++){
            for (int qL = q_min; qL <= q_max; ++qL){
                if (s + qL > win_len) continue;
                if (abs(qL - k) > max_edits) continue;
                AlignScore sc = banded_edit_score(seq + win_start + s, qL, tso_pre, k,
                                                  max_edits, P->n_as_match);
                double er = (k > 0) ? (double)sc.edits / k : 1.0;
                if (sc.edits <= max_edits && er <= P->tso_max_mmr){
                    TSOHit h = { win_start + s, win_start + s + qL,
                                 sc.edits, sc.mm, sc.ins, sc.del, k, er, 1 };
                    th_push(hv, h);
                }
            }
        }
    }
}

/*** FULL: banded edit-distance match with bounded net length shift ***/
static void full_hits(const char *seq, int L, const char *tso, int tL, const TSOParams *P, TSOHits *hv){
    int win_start = (L - P->tail_window > 0) ? L - P->tail_window : 0, win_len = L - win_start;
    int q_min = imax(1, tL - P->tso_max_shift);
    int q_max = tL + P->tso_max_shift;
    for(int i=0; i<win_len; i++){
        for (int qL=q_min; qL<=q_max; ++qL){
            if (i + qL > win_len) continue;
            if (abs(qL - tL) > P->tso_max_mm) continue;
            AlignScore sc = banded_edit_score(seq + win_start + i, qL, tso, tL,
                                              P->tso_max_mm, P->n_as_match);
            if (sc.edits <= P->tso_max_mm){
                int start = win_start + i, end = start + qL;
                TSOHit h={start,end,sc.edits,sc.mm,sc.ins,sc.del,tL,(double)sc.edits/tL,0};
                th_push(hv,h);
            }
        }
    }
}

static void dedup_prune(TSOHits *hv, int min_spacing, int max_keep){
    if (hv->n==0) return; qsort(hv->a, hv->n, sizeof(TSOHit), cmp_desc);
    int w=0; for(int i=0;i<hv->n;i++){
        if(w==0){ hv->a[w++]=hv->a[i]; continue; }
        TSOHit *pr=&hv->a[w-1], *cu=&hv->a[i];
        if (pr->start - cu->start < min_spacing){
            int keep_prev;
            if (pr->partial != cu->partial){
                keep_prev = (pr->partial == 0);
            } else {
                keep_prev = (pr->err_rate < cu->err_rate) ||
                            (pr->err_rate==cu->err_rate && pr->edits <= cu->edits);
            }
            if (!keep_prev) hv->a[w-1] = *cu;
        }else hv->a[w++] = *cu;
    }
    hv->n = w; if (max_keep>0 && hv->n > max_keep) hv->n = max_keep;
}

static inline int better_cut_hit(const TSOHit *a, const TSOHit *b){
    if (!b) return 1;
    if (a->partial != b->partial) return a->partial == 0;
    if (a->start != b->start) return a->start < b->start;
    if (a->err_rate != b->err_rate) return a->err_rate < b->err_rate;
    if (a->edits != b->edits) return a->edits < b->edits;
    return a->end > b->end;
}

static int choose_chain_cut_hit(const TSOHits *hv, int read_len,
                                int max_tail_gap, int max_chain_gap){
    if (hv->n == 0) return -1;

    int best = -1;
    for (int anchor=0; anchor<hv->n; ++anchor){
        int tail_gap = read_len - hv->a[anchor].end;
        if (tail_gap < 0 || tail_gap > max_tail_gap) continue;

        unsigned char *in_chain = (unsigned char*)calloc((size_t)hv->n, 1);
        if (!in_chain) return best;
        in_chain[anchor] = 1;

        int changed = 1;
        while (changed){
            changed = 0;
            for (int i=0; i<hv->n; ++i){
                if (in_chain[i]) continue;
                for (int j=0; j<hv->n; ++j){
                    if (!in_chain[j]) continue;
                    if (hv->a[i].start < hv->a[j].start &&
                        hv->a[i].end + max_chain_gap >= hv->a[j].start){
                        in_chain[i] = 1;
                        changed = 1;
                        break;
                    }
                }
            }
        }

        int chain_pick = anchor;
        for (int i=0; i<hv->n; ++i){
            if (!in_chain[i]) continue;
            if (better_cut_hit(&hv->a[i], &hv->a[chain_pick])) chain_pick = i;
        }

        if (better_cut_hit(&hv->a[chain_pick], best >= 0 ? &hv->a[best] : NULL)){
            best = chain_pick;
        }
        free(in_chain);
    }
    return best;
}

static void tso_scan(const char *seq, int L, const char *tso, int tL, const TSOParams *P, TSOHits *hv){
    th_init(hv); partial_prefix_hits(seq,L,tso,tL,P,hv); full_hits(seq,L,tso,tL,P,hv); dedup_prune(hv,P->min_spacing,P->tso_max_hits);
}

/*************** CLI & options ***************/
typedef struct {
    const char *fastq, *tso, *out_tsv, *out_trim_fastq;
    int tail_window, tso_max_mm, tso_max_shift, tso_min_overlap, tso_max_hits, min_spacing;
    int max_tail_gap, max_chain_gap;
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
"TSOclip " TSOCLIP_VERSION "\n"
"Usage:\n"
"  tsoclip --fastq in.fastq[.gz|-/stdin] --tso SEQ --out-tsv hits.tsv --out-trim-fastq out.fastq[.gz|-/stdout]\n"
"  tsoclip --version\n"
"Scan:\n"
"  --tail-window INT      [150]\n"
"  --tso-max-mm INT       [5]  # max edit distance for full TSO hits\n"
"  --tso-max-shift INT    [4]  # max length shift around TSO/prefix\n"
"  --tso-min-overlap INT  [12]\n"
"  --tso-max-mmr FLOAT    [0.20]\n"
"  --tso-max-hits INT     [100]\n"
"  --min-spacing INT      [6]\n"
"  --max-tail-gap INT     [2]  # max bases after tail-most hit\n"
"  --max-chain-gap INT    [6]  # max gap between chained TSO hits\n"
"Behavior:\n"
"  --n-as-match                 # treat read/TSO 'N' as wildcard [default]\n"
"  --strict-n                   # count 'N' as an ordinary base\n"
"  --concat-prefer-full         # accepted for compatibility; full hits already take priority\n"
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
    o->max_tail_gap=2; o->max_chain_gap=6;
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
        else if(!strcmp(argv[i],"--max-tail-gap") && i+1<argc) o->max_tail_gap=parse_int(argv[++i]);
        else if(!strcmp(argv[i],"--max-chain-gap") && i+1<argc) o->max_chain_gap=parse_int(argv[++i]);
        else if(!strcmp(argv[i],"--threads") && i+1<argc) o->threads=parse_int(argv[++i]);
        else if(!strcmp(argv[i],"--batch-size") && i+1<argc) o->batch_size=parse_int(argv[++i]);
        else if(!strcmp(argv[i],"--concat-prefer-full")) o->concat_prefer_full=1;
        else if(!strcmp(argv[i],"--n-as-match")) o->n_as_match=1;
        else if(!strcmp(argv[i],"--strict-n")) o->n_as_match=0;
        else if(!strcmp(argv[i],"--emit-only-hit") && i+1<argc) o->emit_only_hit=parse_int(argv[++i]);
        else if(!strcmp(argv[i],"--min-keep-len") && i+1<argc) o->min_keep_len=parse_int(argv[++i]);
        else if(!strcmp(argv[i],"--no-json")) o->no_json=1;
        else if(!strcmp(argv[i],"--gzbuf-kb") && i+1<argc) o->gzbuf_kb=parse_int(argv[++i]);
        else if(!strcmp(argv[i],"--gzip-level") && i+1<argc) o->gzip_level=parse_int(argv[++i]);
        else if(!strcmp(argv[i],"--plain-out")) o->plain_out=1;
        else { usage(); return -1; }
    }
    if (!o->fastq || !o->tso || !o->out_tsv || !o->out_trim_fastq){ usage(); return -1; }
    if (o->tail_window <= 0 || o->tso_max_mm < 0 || o->tso_max_shift < 0 ||
        o->tso_min_overlap <= 0 || o->tso_max_mmr < 0.0 || o->batch_size <= 0 ||
        o->min_keep_len < 0 || o->tso_max_hits < 0 || o->min_spacing < 0 ||
        o->max_tail_gap < 0 || o->max_chain_gap < 0 || o->threads < 0 ||
        o->gzbuf_kb < 0 || o->emit_only_hit < 0 || o->emit_only_hit > 1){
        fprintf(stderr, "[ERR] invalid non-positive or negative option\n");
        return -1;
    }
    if (o->tso_max_mmr > 1.0){
        fprintf(stderr, "[ERR] --tso-max-mmr must be <= 1.0\n");
        return -1;
    }
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
    if (argc == 2 && (!strcmp(argv[1], "--version") || !strcmp(argv[1], "-V"))){
        printf("TSOclip %s\n", TSOCLIP_VERSION);
        return 0;
    }
    if (argc == 2 && (!strcmp(argv[1], "--help") || !strcmp(argv[1], "-h"))){
        usage();
        return 0;
    }

    Opt opt; if (parse_args(argc, argv, &opt)!=0) return 1;

    int tL = (int)strlen(opt.tso);
    if (tL <= 0){
        fprintf(stderr, "[ERR] --tso must not be empty\n");
        return 1;
    }
    if (opt.tso_min_overlap >= tL){
        fprintf(stderr, "[ERR] --tso-min-overlap must be smaller than TSO length (%d)\n", tL);
        return 1;
    }
    if (opt.concat_prefer_full){
        fprintf(stderr, "[WARN] --concat-prefer-full is deprecated; full hits already take priority in tail chains\n");
    }
    char *TSO = strdup(opt.tso); for(char *p=TSO; *p; ++p) *p=toupper((unsigned char)*p);
    TSOParams P = { opt.tail_window, opt.tso_max_mm, opt.tso_max_shift,
                    opt.tso_min_overlap, opt.tso_max_hits, opt.min_spacing,
                    opt.max_tail_gap, opt.max_chain_gap,
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

    fprintf(tsv, "read_id\tread_len\ttso_hit_count\tpick_start\tpick_end\tpick_mm\tpick_overlap\tpick_mmr\tpick_partial\tpick_edits\tpick_ins\tpick_del\tpick_error_rate\tcut_start\tcut_end\tcut_len\tcut_reason\tall_hits_json\n");

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
                int pick_idx = choose_chain_cut_hit(&hv, R[i].len, P.max_tail_gap, P.max_chain_gap);
                if (pick_idx >= 0){
                    TSOHit chosen = hv.a[pick_idx];
                    S[i].use_hit=1; S[i].pick=chosen;

                    /* Cut to the earliest retained hit in the chain anchored at
                       the read tail. Full hits take priority over partial hits
                       when both are present in that chain. */
                    S[i].out_len = chosen.start;
                } else {
                    memset(&S[i].pick, 0, sizeof(S[i].pick));
                }

                if (!opt.no_json){
                    int est = hv.n*96 + 16; S[i].all_json = (char*)malloc(est);
                    char *w = S[i].all_json; int rem=est;
                    w += snprintf(w, rem, "["); rem = est - (w - S[i].all_json);
                    for(int k=0;k<hv.n;k++){
                        w += snprintf(w, rem, "%s{\"s\":%d,\"e\":%d,\"ed\":%d,\"mm\":%d,\"ins\":%d,\"del\":%d,\"ov\":%d,\"er\":%.4f,\"p\":%d}",
                                      (k?",":""), hv.a[k].start, hv.a[k].end,
                                      hv.a[k].edits, hv.a[k].mm, hv.a[k].ins, hv.a[k].del,
                                      hv.a[k].overlap, hv.a[k].err_rate, hv.a[k].partial);
                        rem = est - (w - S[i].all_json);
                        if (rem < 128){ int used=(int)(w - S[i].all_json); est<<=1; S[i].all_json=(char*)realloc(S[i].all_json, est); w=S[i].all_json+used; rem=est-used; }
                    }
                    snprintf(w, rem, "]");
                } else {
                    S[i].all_json = strdup("[]");
                }
            } else {
                memset(&S[i].pick, 0, sizeof(S[i].pick));
                S[i].all_json = opt.no_json ? strdup("[]") : strdup("[]");
            }
            th_free(&hv);
        }

        for(int i=0;i<n;i++){
            if (S[i].hit_count>0 && S[i].use_hit){
                fprintf(tsv, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.3f\t%s\t%d\t%d\t%d\t%.3f\t%d\t%d\t%d\t%s\t%s\n",
                    R[i].id, R[i].len, S[i].hit_count,
                    S[i].pick.start, S[i].pick.end, S[i].pick.mm, S[i].pick.overlap,
                    S[i].pick.err_rate, S[i].pick.partial? "YES":"NO",
                    S[i].pick.edits, S[i].pick.ins, S[i].pick.del, S[i].pick.err_rate,
                    S[i].out_len, R[i].len, R[i].len - S[i].out_len,
                    S[i].pick.partial? "PARTIAL":"FULL", S[i].all_json);
            } else if (S[i].hit_count>0 && !S[i].use_hit){
                fprintf(tsv, "%s\t%d\t%d\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tNO_CUT\t%s\n",
                    R[i].id, R[i].len, S[i].hit_count, S[i].all_json);
            } else {
                fprintf(tsv, "%s\t%d\t0\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tNO_HIT\t[]\n", R[i].id, R[i].len);
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
