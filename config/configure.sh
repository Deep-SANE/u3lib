#!/usr/bin/env bash
# u3lib configuration — resolve every build dependency and write
# config/config.generated.mk. Run via `make config` (or directly).
#
# Each dependency is resolved by the chain:
#   1. detect  — standard locations (+ pkg-config / Homebrew where applicable)
#   2. ask      — if not detected, prompt for the path to an existing install
#   3. install  — if none given, download the pinned upstream source + build it
#                 into external-libs/ (which is not part of the repository)
# The compiler and OpenMP are required system tools and are detected, not installed.
#
# Non-interactive use (CI, scripted installs):
#   ./configure.sh --yes                  detect, and install anything undetected
#   ./configure.sh --cxx=g++-16 --boost=/path --eigen=/path/eigen3 \
#                  --fmt=/path --wignernj=fetch --krans=/path/to/krans-install --yes
# A flag value of "fetch" forces step 3 (download + build) for that dependency;
# a path forces step 2 (use that existing install).
set -uo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
EXT="$ROOT/external-libs"
GEN="$ROOT/config/config.generated.mk"
JOBS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)"

# Option overrides are stored in plain variables OPT_<key> (no associative
# arrays, for bash 3.2 compatibility on stock macOS).
YES=0
OPT_cxx=""; OPT_cc=""; OPT_boost=""; OPT_eigen=""; OPT_fmt=""
OPT_wignernj=""; OPT_krans=""
for arg in "$@"; do
  case "$arg" in
    --yes|-y) YES=1 ;;
    --jobs=*) JOBS="${arg#*=}" ;;
    --cxx=*)  OPT_cxx="${arg#*=}" ;;
    --cc=*)   OPT_cc="${arg#*=}" ;;
    --boost=*) OPT_boost="${arg#*=}" ;;
    --eigen=*) OPT_eigen="${arg#*=}" ;;
    --fmt=*)  OPT_fmt="${arg#*=}" ;;
    --wignernj=*) OPT_wignernj="${arg#*=}" ;;
    --krans=*)  OPT_krans="${arg#*=}" ;;
    --help|-h)
      sed -n '2,20p' "$0"; exit 0 ;;
    *) echo "unknown option: $arg" >&2; exit 2 ;;
  esac
done

say()  { printf '\033[1;34m==>\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33mnote:\033[0m %s\n' "$*" >&2; }
die()  { printf '\033[1;31merror:\033[0m %s\n' "$*" >&2; exit 1; }

# ask KEY "prompt" "default" -> prints the chosen value.
# Honors --KEY=... overrides and --yes (take the default without prompting).
ask() {
  local key="$1" prompt="$2" def="$3" ans ov
  eval "ov=\${OPT_$key:-}"
  if [[ -n "$ov" ]]; then echo "$ov"; return; fi
  if [[ $YES -eq 1 ]]; then echo "$def"; return; fi
  read -r -p "$prompt [$def]: " ans </dev/tty || ans=""
  echo "${ans:-$def}"
}

first_in_path() { local c; for c in "$@"; do if command -v "$c" >/dev/null 2>&1; then command -v "$c"; return; fi; done; }
brew_prefix()   { command -v brew >/dev/null 2>&1 && brew --prefix "$1" 2>/dev/null; }

# download URL OUTFILE -- cross-platform (curl or wget).
download() {
  if command -v curl >/dev/null 2>&1; then curl -fsSL "$1" -o "$2"
  elif command -v wget >/dev/null 2>&1; then wget -qO "$2" "$1"
  else die "need curl or wget to download $1"; fi
}

have_pkgconfig() { command -v pkg-config >/dev/null 2>&1; }
is_install()     { [[ "$1" == INSTALL || "$1" == fetch || "$1" == build ]]; }

# Resolve one dependency following the chain:
#   --KEY=<path|fetch> flag  >  auto-detected  >  ask for a path  >  install.
# Prints the chosen path, or the literal "INSTALL" when it should be built.
# (Blank answer, or --yes with nothing detected, means install.)
decide() {
  local key="$1" name="$2" detected="$3" ov ans
  eval "ov=\${OPT_$key:-}"
  if [[ -n "$ov" ]]; then echo "$ov"; return; fi
  if [[ -n "$detected" ]]; then say "found $name: $detected" >&2; echo "$detected"; return; fi
  if [[ $YES -eq 1 ]]; then echo "INSTALL"; return; fi
  read -r -p "$name not detected. Path to an existing install, or Enter to download & build: " ans </dev/tty || ans=""
  [[ -z "$ans" ]] && echo "INSTALL" || echo "$ans"
}

# Detection probes: standard locations, plus pkg-config where it applies
# (covers many Linux distro and HPC/module setups). Echo a valid path or nothing.
detect_boost() {
  local p
  for p in "$BREW_BOOST" /opt/homebrew/opt/boost /usr/local /usr; do
    [[ -n "$p" && -d "$p/include/boost" ]] && { echo "$p"; return; }
  done
}
detect_eigen() {
  local inc p
  if have_pkgconfig && pkg-config --exists eigen3 2>/dev/null; then
    inc="$(pkg-config --variable=includedir eigen3 2>/dev/null)"
    [[ -n "$inc" && -f "$inc/eigen3/Eigen/Dense" ]] && { echo "$inc/eigen3"; return; }
    [[ -n "$inc" && -f "$inc/Eigen/Dense" ]] && { echo "$inc"; return; }
  fi
  for p in "${BREW_EIGEN:+$BREW_EIGEN/include/eigen3}" \
           /opt/homebrew/opt/eigen/include/eigen3 /usr/include/eigen3 /usr/local/include/eigen3; do
    [[ -n "$p" && -f "$p/Eigen/Dense" ]] && { echo "$p"; return; }
  done
}
detect_fmt() {
  local inc p
  if have_pkgconfig && pkg-config --exists fmt 2>/dev/null; then
    inc="$(pkg-config --variable=includedir fmt 2>/dev/null)"
    [[ -n "$inc" && -f "$inc/fmt/core.h" ]] && { echo "${inc%/include}"; return; }
  fi
  for p in "$BREW_FMT" /opt/homebrew/opt/fmt /usr/local /usr; do
    [[ -n "$p" && -f "$p/include/fmt/core.h" ]] && { echo "$p"; return; }
  done
}

# ---------------------------------------------------------------------------
say "Configuring u3lib (root: $ROOT)"

# --- C++ compiler (needs C++20; real GCC or clang with OpenMP) ---
CXX_DEF="$(first_in_path g++-16 g++-15 g++-14 g++ c++ || true)"; CXX_DEF="${CXX_DEF:-g++}"
CXX="$(ask cxx "C++ compiler (C++20)" "$CXX_DEF")"
command -v "$CXX" >/dev/null 2>&1 || die "C++ compiler not found: $CXX"

# --- C compiler (only used if libwignernj is built locally) ---
CC_DEF="$(first_in_path gcc-16 gcc-15 gcc-14 gcc cc || true)"; CC_DEF="${CC_DEF:-cc}"
CC="$(ask cc "C compiler (for a local libwignernj build)" "$CC_DEF")"

# Precompute Homebrew prefixes once (empty if brew / the formula is absent).
BREW_BOOST="$(brew_prefix boost || true)"
BREW_EIGEN="$(brew_prefix eigen || true)"
BREW_FMT="$(brew_prefix fmt || true)"

# --- Boost (headers; BOOST_ROOT must contain include/boost) ---
BOOST_VER=1.86.0; BOOST_USCORE=1_86_0
BOOST_ROOT="$(decide boost "Boost (headers)" "$(detect_boost)")"
if is_install "$BOOST_ROOT"; then
  say "Downloading Boost $BOOST_VER headers -> external-libs/boost (large, ~130 MB)"
  command -v tar >/dev/null 2>&1 || die "tar required to unpack Boost"
  rm -rf "$EXT/boost"; mkdir -p "$EXT/boost/include"
  download "https://archives.boost.io/release/$BOOST_VER/source/boost_$BOOST_USCORE.tar.gz" "$EXT/boost/boost.tar.gz" \
    || die "Boost download failed"
  tar -xzf "$EXT/boost/boost.tar.gz" -C "$EXT/boost" "boost_$BOOST_USCORE/boost" || die "Boost extract failed"
  mv "$EXT/boost/boost_$BOOST_USCORE/boost" "$EXT/boost/include/boost"
  rm -rf "$EXT/boost/boost_$BOOST_USCORE" "$EXT/boost/boost.tar.gz"
  BOOST_ROOT="$EXT/boost"
fi
[[ -d "$BOOST_ROOT/include/boost" ]] || die "Boost headers not found under $BOOST_ROOT/include/boost"

# --- Eigen (header-only; EIGEN_DIR contains Eigen/) ---
EIGEN_DIR="$(decide eigen "Eigen3" "$(detect_eigen)")"
if is_install "$EIGEN_DIR"; then
  say "Fetching Eigen 3.4.0 (header-only) -> external-libs/eigen"
  command -v git >/dev/null 2>&1 || die "git required to fetch Eigen"
  rm -rf "$EXT/eigen"
  git clone --depth 1 --branch 3.4.0 https://gitlab.com/libeigen/eigen.git "$EXT/eigen" >/dev/null 2>&1 \
    || die "Eigen fetch failed"
  EIGEN_DIR="$EXT/eigen"
fi
[[ -f "$EIGEN_DIR/Eigen/Dense" ]] || die "Eigen not found under $EIGEN_DIR (expected Eigen/Dense)"

# --- fmt (header-only; FMT_DIR contains include/fmt) ---
FMT_TAG=10.2.1
FMT_DIR="$(decide fmt "fmt" "$(detect_fmt)")"
if is_install "$FMT_DIR"; then
  say "Fetching fmt $FMT_TAG (header-only) -> external-libs/fmt"
  command -v git >/dev/null 2>&1 || die "git required to fetch fmt"
  rm -rf "$EXT/fmt"
  git clone --depth 1 --branch "$FMT_TAG" https://github.com/fmtlib/fmt "$EXT/fmt" >/dev/null 2>&1 \
    || die "fmt clone failed"
  FMT_DIR="$EXT/fmt"
fi
[[ -f "$FMT_DIR/include/fmt/core.h" ]] || die "fmt not found under $FMT_DIR (expected include/fmt/core.h)"

# Pinned upstream source for the vendored-at-build-time numerical library.
WIGNERNJ_REPO="https://github.com/susilehtola/libwignernj"
WIGNERNJ_TAG="v0.7.0"

# --- libwignernj (existing install dir, or 'fetch' to download + build) ---
valid_wignernj() { [[ -f "$1/include/wignernj.hpp" ]] && ls "$1"/lib/libwignernj.* >/dev/null 2>&1; }
WJ_DET=""; valid_wignernj "$EXT/libwignernj/install" && WJ_DET="$EXT/libwignernj/install"
WIGNERNJ_DIR="$(decide wignernj "libwignernj" "$WJ_DET")"
if is_install "$WIGNERNJ_DIR"; then
  say "Fetching libwignernj $WIGNERNJ_TAG and building (C only) -> external-libs/libwignernj/install"
  command -v git >/dev/null 2>&1 || die "git required to fetch libwignernj"
  command -v cmake >/dev/null 2>&1 || die "cmake required to build libwignernj"
  rm -rf "$EXT/libwignernj"
  git clone --depth 1 --branch "$WIGNERNJ_TAG" "$WIGNERNJ_REPO" "$EXT/libwignernj/src" >/dev/null 2>&1 \
    || die "libwignernj clone failed ($WIGNERNJ_REPO @ $WIGNERNJ_TAG)"
  cmake -S "$EXT/libwignernj/src" -B "$EXT/libwignernj/build" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER="$CC" \
    -DBUILD_FORTRAN=OFF -DBUILD_TESTS=OFF -DBUILD_EXAMPLES=OFF \
    -DBUILD_CXX_TESTS=OFF -DBUILD_PYTHON=OFF -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_INSTALL_PREFIX="$EXT/libwignernj/install" >/dev/null || die "libwignernj: cmake configure failed"
  cmake --build "$EXT/libwignernj/build" -j "$JOBS" >/dev/null || die "libwignernj: build failed"
  cmake --install "$EXT/libwignernj/build" >/dev/null || die "libwignernj: install failed"
  WIGNERNJ_DIR="$EXT/libwignernj/install"
fi
valid_wignernj "$WIGNERNJ_DIR" || die "libwignernj not usable at $WIGNERNJ_DIR (need include/wignernj.hpp + lib/libwignernj.*)"

# --- krans (Deep-SANE/krans): sibling checkout; auto-cloned + built here if missing ---
KRANS_REPO="git@github.com:Deep-SANE/krans.git"  # TODO: switch to public https URL once KraNS is published
valid_krans() { [[ -d "$1/include/krans" && -f "$1/build/libkrans_core.a" && -f "$1/build/libkrans_eigen.a" ]]; }
# Location: --krans=PATH, env KRANS_DIR, or the sibling ../krans next to u3lib.
KRANS_DIR="${OPT_krans:-${KRANS_DIR:-}}"
[[ -z "$KRANS_DIR" ]] && KRANS_DIR="$(cd "$ROOT/.." && pwd)/krans"
# Clone it next to u3lib if the source is not present yet (helps users who
# forgot to install it): a plain clone, not vendored under external-libs.
if [[ ! -d "$KRANS_DIR/include/krans" ]]; then
  command -v git >/dev/null 2>&1 || die "git required to fetch krans"
  say "krans not found; cloning -> $KRANS_DIR"
  git clone "$KRANS_REPO" "$KRANS_DIR" >/dev/null 2>&1 || die "krans clone failed ($KRANS_REPO).
Check your SSH access to the Deep-SANE organization, or clone it manually:
  git clone $KRANS_REPO \"$KRANS_DIR\""
fi
# Build the static libraries in place if they are not present yet.
if [[ ! -f "$KRANS_DIR/build/libkrans_core.a" || ! -f "$KRANS_DIR/build/libkrans_eigen.a" ]]; then
  command -v cmake >/dev/null 2>&1 || die "cmake required to build krans"
  say "Building krans in place -> $KRANS_DIR/build"
  cmake -S "$KRANS_DIR" -B "$KRANS_DIR/build" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER="$CXX" >/dev/null || die "krans: cmake configure failed"
  cmake --build "$KRANS_DIR/build" -j >/dev/null || die "krans: build failed"
fi
valid_krans "$KRANS_DIR" || die "krans not usable at $KRANS_DIR (need include/krans + build/libkrans_{core,eigen}.a)"
say "Found krans at $KRANS_DIR"

# --- OpenMP (probe the chosen compiler) ---
USE_OPENMP=0
if printf 'int main(){return 0;}\n' | "$CXX" -fopenmp -x c++ - -o /dev/null >/dev/null 2>&1; then
  USE_OPENMP=1
else
  warn "compiler $CXX does not accept -fopenmp; building without OpenMP"
fi

# ---------------------------------------------------------------------------
say "Writing $GEN"
{
  echo "# Generated by config/configure.sh on $(date). Re-run \`make config\` to change."
  echo "CXX := $CXX"
  echo "CXXFLAGS := -std=c++20 -O3"
  echo "LDFLAGS :="
  echo "USE_WIGNERNJ := 1"
  echo "WIGNERNJ_DIR := $WIGNERNJ_DIR"
  echo "USE_KRANS := 1"
  echo "KRANS_DIR := $KRANS_DIR"
  echo "NULLSPACE_SOLVER := krans"
  echo "BOOST_ROOT := $BOOST_ROOT"
  echo "EIGEN_DIR := $EIGEN_DIR"
  if [[ "$FMT_DIR" != "none" ]]; then
    echo "HAVE_FMT := 1"
    echo "FMT_HEADER_ONLY := 1"
    echo "FMT_DIR := $FMT_DIR"
  fi
  if [[ "$USE_OPENMP" -eq 1 ]]; then
    echo "USE_OPENMP := 1"
    echo "OMP_CXXFLAGS := -fopenmp"
    echo "OMP_LDFLAGS := -fopenmp"
  fi
} > "$GEN"

say "Done. Configuration summary:"
cat "$GEN"
echo
say "Next: make  (then: make test, make python)"
