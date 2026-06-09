# Release Checklist

1. Update `VERSION` and `CHANGELOG.md`.
2. Run `make clean && make && make check`.
3. Build the container with `docker build -t tsoclip:$(cat VERSION) .`.
4. Create a tagged GitHub release from a clean commit.
5. Archive the release on Zenodo and record the DOI in the README or manuscript.
6. Re-run manuscript examples against the archived release, not a moving branch.
