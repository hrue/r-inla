# The binaries, what is inside them, and what pairs with what

Versions that float (the compiler, ARMPL, OpenBLAS) are set in
[`toolchain.env`](toolchain.env) and recorded per-artifact in the `BUILDINFO`
file each archive carries. This table is the map; `BUILDINFO` is the ground
truth for any given artifact.

## inla bundles (GitHub Releases on this repository)

| bundle | runs on | CPU target | compiler | BLAS/LAPACK | libstiles inside | notes |
|---|---|---|---|---|---|---|
| `inla-linux-x86_64-portable` | any Linux, glibc 2.28+ (RHEL 8+), **AVX2 CPU** | baseline x86-64 + runtime dispatch (CLONE_TARGETS, MKL) | newest gcc-toolset (manylinux) | Intel MKL, static, GNU-threaded (shares libgomp) | `libstiles-linux-x86_64` | the AVX2 floor comes from libstiles, not inla |
| `inla-linux-x86_64-v3-portable` | Linux, glibc 2.38+, 2015+ CPU | `-march=x86-64-v3` everywhere | newest GCC (PPA) | Intel MKL, static | `libstiles-linux-x86_64-v3-mkl` | the performance build |
| `inla-linux-arm64-armv82-portable` | Linux arm64, glibc 2.38+ (Graviton2+, Ampere, RPi 5) | `-march=armv8.2-a -mtune=neoverse-n1` | newest GCC (PPA) | ARMPL, static, serial | `libstiles-linux-arm64-armv82-armpl` | the only arm64 bundle: the baseline-ISA arm builds hit a latent heap bug (see repo history) |
| `inla-macos-arm64-portable` | macOS (floor set by Homebrew libs; see BUILDINFO) | `-mcpu=apple-m1` | Homebrew GCC | ARMPL, serial, bundled | `libstiles-macos-arm64-gcc-armpl` | rgeneric loads the running machine's R |
| `inla-macos-x86_64-portable` | macOS Intel (floor per BUILDINFO) | `-mtune=generic` | Homebrew GCC | OpenBLAS, serial locking-safe, bundled | `libstiles-macos-intel-x86_64-gcc-openblas` | Accelerate rejected: its threading cannot be controlled |
| `inla-windows-x86_64.zip` | Windows 10+, with or without R | `-mtune=generic` | MinGW-w64 GCC | R's Rblas/Rlapack, bundled (whole R DLL closure) | `libstiles-windows-x86_64` | mimalloc allocator; libstiles.dll linked directly (its implib was defective before 2026.8.19) |

Common to every bundle: `WITH_LIBR=2` (no libR linked; rgeneric dlopens the
running machine's R), one OpenMP runtime (GNU libgomp), static libstdc++,
LTO, external model packages compiled in, `BUILDINFO` at the archive root.

## libstiles artifacts (GitHub Releases on esmail-abdulfattah/sTiles)

The first five names are fixed forever: pip and R download them by exact name.

| artifact | toolchain | CPU target | BLAS/LAPACK (embedded unless noted) | OpenMP | consumer |
|---|---|---|---|---|---|
| `libstiles-linux-x86_64` | gcc-toolset (manylinux, glibc 2.27) | AVX2 (haswell) | MKL, sequential | libgomp | pip/R, inla classic |
| `libstiles-linux-arm64` | gcc-toolset (manylinux) | armv8 baseline | OpenBLAS, serial | libgomp | pip/R |
| `libstiles-macos-apple-arm64` | clang | apple-m1 | Accelerate (system) | **libomp** (LLVM) | pip/R only — never link into a GCC build |
| `libstiles-macos-intel-x86_64` | clang | generic | Accelerate (system) | **libomp** (LLVM) | pip/R only |
| `libstiles-windows-x86_64` | MSYS2 UCRT GCC | AVX2 | OpenBLAS, serial (DLL, bundled) | libgomp | pip/R, inla |
| `libstiles-linux-x86_64-v3-mkl` | newest GCC | AVX2/x86-64-v3 | MKL | libgomp | inla v3 bundle |
| `libstiles-linux-arm64-armv8-armpl` | gcc-toolset | armv8 baseline | ARMPL, serial | libgomp | spare (was inla classic arm) |
| `libstiles-linux-arm64-armv82-armpl` | newest GCC | armv8.2 | ARMPL, serial | libgomp | inla arm64 bundle |
| `libstiles-macos-arm64-gcc-armpl` | Homebrew GCC | apple-m1 | ARMPL (external, from the inla bundle) | libgomp | inla mac arm |
| `libstiles-macos-intel-x86_64-gcc-openblas` | Homebrew GCC | generic | OpenBLAS, serial (bundled + GCC runtime) | libgomp | inla mac Intel |
| `libstiles-macos-*-static` | clang | native | none (bring your own CBLAS+LAPACKE) | none | integrators |

## The pairing rule

An inla bundle and its libstiles must agree on the **OpenMP runtime** (both
GNU libgomp — the clang/libomp pip assets must never be linked into a GCC
inla), and the libstiles must not be built by a **newer** compiler than the
bundle's runtime libraries. Everything else (BLAS choice, ISA level) is
per-artifact and encoded in the names above. Both projects read
`toolchain.env` from this repository, so the versions cannot drift silently;
compare the `BUILDINFO` files of a bundle and its libstiles to verify a pair.
