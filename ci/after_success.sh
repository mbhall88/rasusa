# This script takes care of uploading code coverage to codecov

set -ex

main() {
  wget https://github.com/SimonKagstrom/kcov/archive/master.tar.gz &&
  tar xzf master.tar.gz &&
  cd kcov-master &&
  mkdir build &&
  cd build &&
  cmake .. &&
  make &&
  make install DESTDIR=../../kcov-build &&
  cd ../.. &&
  rm -rf kcov-master
  export PATH=$(realpath ./kcov-build/usr/local/bin):"$PATH"
  for file in target/debug/"$PROJECT_NAME"-*; do [ -x "${file}" ] || continue; mkdir -p "target/cov/$(basename $file)"; kcov --exclude-pattern=/.cargo,/usr/lib --verify "target/cov/$(basename $file)" "$file" || exit 1; done &&
  bash <(curl -s https://codecov.io/bash) &&
  echo "Uploaded code coverage"
}

# we don't run kcov on osx and only run on stable
if [[ -z "$TRAVIS_OS_NAME" = "linux" && "$TRAVIS_RUST_VERSION" = "stable" ]]; then
    main
fi
