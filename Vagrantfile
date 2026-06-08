# -*- mode: ruby -*-
# vi: set ft=ruby :

# Vagrantfile for testing Rpcop with ASAN/UBSAN on Linux

Vagrant.configure("2") do |config|
  # Use bento/ubuntu-22.04, matching the host architecture when needed.
  config.vm.box = "bento/ubuntu-22.04"
  config.vm.hostname = "rpcop-sanitizer"

  # Use rsync for shared folders (doesn't require interactive input)
  config.vm.synced_folder ".", "/vagrant", type: "rsync", rsync__exclude: [".git/", "*.o", "*.so", "*.tar.gz"]

  # VirtualBox provider
  config.vm.provider "virtualbox" do |vb|
    vb.memory = "4096"
    vb.cpus = 2
  end

  # Provisioning script - install R and dependencies, then run sanitizer tests
  config.vm.provision "shell", inline: <<-SHELL
    set -eu

    echo ">>> Updating apt repositories..."
    apt-get update

    echo ">>> Installing build dependencies..."
    apt-get install -y \
      build-essential \
      gfortran \
      libcurl4-openssl-dev \
      libssl-dev \
      libxml2-dev \
      libuv1-dev \
      libfontconfig1-dev \
      libharfbuzz-dev \
      libfribidi-dev \
      libfreetype6-dev \
      libpng-dev \
      libtiff5-dev \
      libjpeg-dev \
      wget \
      ca-certificates \
      curl \
      libasan6 \
      libubsan1 \
      libasan8

    echo ">>> Adding CRAN repository for latest R..."
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
    UBUNTU_ARCH=$(dpkg --print-architecture)
    echo "deb [arch=$UBUNTU_ARCH] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" | tee /etc/apt/sources.list.d/cran-r.list
    apt-get update

    echo ">>> Installing R..."
    apt-get install -y r-base r-base-dev

    cd /vagrant

    echo ">>> Installing R package dependencies..."
    mkdir -p /tmp/Rlibs
    R_PROFILE_USER=/dev/null R_LIBS=/tmp/Rlibs Rscript tools/install-sanitizer-deps.R

    echo ""
    echo ">>> Building Rpcop with ASAN/UBSAN..."
    rm -f src/*.o src/*.so

    # Build with sanitizers enabled
    # Use --no-test-load to avoid ASAN preload issues during install
    R_PROFILE_USER=/dev/null R_LIBS=/tmp/Rlibs USE_SANITIZER=1 R CMD INSTALL --preclean --no-multiarch --no-test-load --library=/tmp/Rlibs . 2>&1

    echo ""
    echo ">>> Finding ASAN library..."
    ASAN_LIB=$(find /usr -name "libasan.so*" 2>/dev/null | head -1)
    echo "ASAN library: $ASAN_LIB"
    if [ -z "$ASAN_LIB" ]; then
      echo "Could not find libasan.so."
      exit 1
    fi

    echo ""
    echo ">>> Running tests with sanitizers..."
    export LD_PRELOAD="$ASAN_LIB"
    export ASAN_OPTIONS="detect_leaks=0:detect_stack_use_after_return=1:halt_on_error=1:detect_container_overflow=0"
    export UBSAN_OPTIONS="print_stacktrace=1:halt_on_error=1"
    if ! R_PROFILE_USER=/dev/null R_LIBS=/tmp/Rlibs /usr/lib/R/bin/Rscript tools/run-sanitizer-tests.R > /tmp/sanitizer_output.txt 2>&1; then
      cat /tmp/sanitizer_output.txt
      exit 1
    fi
    cat /tmp/sanitizer_output.txt

    echo ""
    echo ">>> Checking for sanitizer errors in Rpcop code..."
    if grep -E "ERROR: (Address|Leak)Sanitizer|runtime error:" /tmp/sanitizer_output.txt; then
      echo "Sanitizer errors detected."
      exit 1
    fi
    echo "No sanitizer errors detected."
  SHELL
end
