# Checking the package on the platforms this machine is not

The compiled kernels are written to return values identical bit for bit to the
R implementations they replace. Two things that identity depends on cannot be
observed on an Apple Silicon Mac:

* **The accumulator.** R sums in `LDOUBLE`, which is `long double` where
  `capabilities("long.double")` is TRUE. It is FALSE here, and on aarch64
  `long double` is `double`, so the wide branch of every reduction template is
  compiled and never selected.
* **The BLAS.** Everything here goes to Accelerate.

Three containers cover what is missing. Each was found to matter: the first two
each exposed a defect that no run on this machine could have.

| file | what it is | why |
| --- | --- | --- |
| `Dockerfile.linux` | Debian, aarch64 | `long double` is 128-bit and the BLAS is the reference one |
| `Dockerfile.asan`  | CRAN's `clang-asan` flavour, x86_64 | R itself built with `-fsanitize=address,undefined`; also 80-bit `long double` |
| `Dockerfile.nold`  | CRAN's no-long-double flavour, x86_64 | the platform has a wide `long double` and R is built not to use it, which no compile-time test can tell from a platform that has none |

AddressSanitizer cannot be run on macOS: its runtime has to be present before
the process starts, and injecting it with `DYLD_INSERT_LIBRARIES` hangs in
dyld's own initialiser before `main`.

## Running them

```sh
colima start --vm-type=vz --vz-rosetta --cpu 6 --memory 10 --disk 80

docker build -f Dockerfile.linux -t drbayes-linux .
docker build --platform linux/amd64 -f Dockerfile.asan -t drbayes-asan .
docker build --platform linux/amd64 -f Dockerfile.nold -t drbayes-nold .

# then, for each image, with the package mounted at /pkg:
#   R CMD INSTALL /pkg && NOT_CRAN=true R_PROFILE_USER=/dev/null Rscript run-suite.R
```

`NOT_CRAN=true` matters: the tests that compare the compiled kernels with the R
they replace carry `skip_on_cran()`, and they are the point of the exercise.
`R_PROFILE_USER=/dev/null` matters too, or `renv/activate.R` resets
`.libPaths()` and the run silently uses a different build of the package.
