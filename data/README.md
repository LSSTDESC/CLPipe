# CLPipe/data

This folder is intentionally empty in git — only this file is tracked. The
actual chains, sacc files, and mock realizations referenced by the example
notebooks (as `../../data/<example>/...`) are hosted externally on the
NERSC Science Gateway portal, not committed to the repo.

**Download instructions for users: see the main [README.md](../README.md).**

## Where it lives

- Portal (read-only, no auth): https://portal.nersc.gov/cfs/lsst/clpipe/data/
- Filesystem (NERSC CFS): `/global/cfs/cdirs/lsst/www/clpipe/data/`

## Structure

```
data/
  cosmodc2_redmapper/
    chains/    *.fits snapshots (see clpipe.utils.cosmosis_mcmc_plots.chains_to_fits)
    sacc/      *.sacc covariance files
  capish_simulation/
    chains/
    *.sacc
    data_generation/mocks_seeds/   ~500MB of mock seed realizations
```

Keep this layout in sync with the `../../data/...` paths hardcoded in the
example notebooks/scripts if you ever reorganize it.

## Refreshing the hosted bundle from CCIN2P3

If files are stored in CC-IN2P3 or another CC (`/sps/lsst/users/ebarroso/CLPipe/examples/`).
After regenerating chains/sacc there (the notebooks already save straight to
`../../data/...` locally):

1. **Stage a clean copy** on CC-IN2P3, mirroring this structure. Hardlinks
   are instant and use no extra disk space (same filesystem):
   ```bash
   cp -al examples/cosmodc2_redmapper/data/... /path/to/staging/...
   ```
   (or just point the transfer directly at the local `data/` tree that the
   notebooks already populate, if it only contains what should ship)

2. **Transfer to NERSC.** Globus is recommended for CC-IN2P3 -> NERSC
   facility-to-facility transfers (one-time browser auth, then
   resumable/checksummed). Direct `rsync` also works if you set up an SSH
   key on the NERSC side, added to CC-IN2P3's `~/.ssh/authorized_keys`.

3. **Fix permissions.** `rsync -a` preserves CC-IN2P3's restrictive
   source permissions (group-only), which breaks the portal (403 /
   "unable to read htaccess file"). Every sync, re-open the tree:
   ```bash
   chmod -R a+rX /global/cfs/cdirs/lsst/www/clpipe/data
   ```

4. Verify: https://portal.nersc.gov/cfs/lsst/clpipe/data/ should list both
   subdirectories.
