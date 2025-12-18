#!/bin/bash
set -e
export LC_ALL="C" # suppress warnings in p4est installation

# Set default values
INTERACTIVE=true
SKIP_PREREQUISITE=false
PREREQUISITES="gcc make cmake open-mpi hdf5-mpi boost zlib lapack tbb (assimp)"
PREREQUISITES_BREW="gcc make cmake open-mpi hdf5-mpi boost lapack tbb assimp"
PREREQUISITES_APT="build-essential gcc make cmake openmpi-bin libhdf5-openmpi-dev libboost-all-dev  zlib1g-dev liblapack-dev libtbb2 libtbb2-dev libassimp-dev"
if [[ $(uname) == "Darwin" ]]; then
    PACKAGE_INSTALLER="brew"
else
    PACKAGE_INSTALLER="apt"
fi

INSTALL_ASSIMP=false
ASSIMP_DIR=${ASSIMP_DIR:-"$HOME/.local/lib/assimp"}
# Check here for new versions: https://github.com/assimp/assimp/releases
ASSIMP_VERSION="6.0.2"

INSTALL_PETSC=false
PETSC_DIR=${PETSC_DIR:-"$HOME/.local/lib/petsc"}
# Check here for new versions: https://gitlab.com/petsc/petsc/-/tags
PETSC_VERSION="3.24.2"
if [[ $(uname) == "Darwin" ]]; then
    PETSC_ARCH=${PETSC_ARCH:-"arch-darwin-c-debug"}
else
    PETSC_ARCH=${PETSC_ARCH:-"arch-linux-c-debug"}
fi
PETSC_FLAGS=${PETSC_FLAGS:-"--download-scalapack --download-mumps --download-hypre"}

INSTALL_P4EST=false
P4EST_DIR=${P4EST_DIR:-"$HOME/.local/lib/p4est"}
# Check here for new versions: https://github.com/cburstedde/p4est/tags
P4EST_VERSION="2.8.7"

INSTALL_DEAL_II=false
DEAL_II_DIR=${DEAL_II_DIR:-"$HOME/.local/lib/dealii"}
# Check here for new versions: https://github.com/dealii/dealii/releases
DEAL_II_VERSION="9.7.1"
DEAL_II_FLAGS=${DEAL_II_FLAGS:-""}

if [[ $(uname) == "Darwin" ]]; then
    NUMBER_JOBS=$(sysctl -n hw.ncpu)
else
    NUMBER_JOBS=$(nproc)
fi
RUN_TEST=false
WORKPWD=${WORKPWD:-$(pwd)}

# Define functions
function set_configuration {
    # Prompt the user for input
    read -rp "Is assimp already installed on your machine? [y/N]: " install_assimp_input
    read -rp "Enter the installation directory for assimp (default: $ASSIMP_DIR): " assimp_dir_input
    read -rp "Is PETSc already installed on your machine? [y/N]: " install_petsc_input
    read -rp "Enter the installation directory for PETSc (default: $PETSC_DIR): " petsc_dir_input
    read -rp "Is p4est already installed on your machine? [y/N]: " install_p4est_input
    read -rp "Enter the installation directory for p4est (default: $P4EST_DIR): " p4est_dir_input
    read -rp "Is deal.II already installed on your machine? [y/N]: " install_deal_ii_input
    read -rp "Enter the installation directory for deal.II (default: $DEAL_II_DIR): " deal_ii_dir_input
    read -rp "Enter the number of processors to use for compilation (default: $NUMBER_JOBS): " number_jobs_input
    read -rp "Do you want to run the test suite after installation? [y/N]: " run_test_input

    # Set the configuration based on user input
    if [[ $install_assimp_input =~ ^[Nn]$ ]]; then
        INSTALL_ASSIMP=true
    fi
    if [[ $assimp_dir_input ]]; then
        ASSIMP_DIR=$assimp_dir_input
    fi
    if [[ $install_petsc_input =~ ^[Nn]$ ]]; then
        INSTALL_PETSC=true
    fi
    if [[ $petsc_dir_input ]]; then
        PETSC_DIR=$petsc_dir_input
    fi
    if [[ $install_p4est_input =~ ^[Nn]$ ]]; then
        INSTALL_P4EST=true
    fi
    if [[ $p4est_dir_input ]]; then
        P4EST_DIR=$p4est_dir_input
    fi
    if [[ $install_deal_ii_input =~ ^[Nn]$ ]]; then
        INSTALL_DEAL_II=true
    fi
    if [[ $deal_ii_dir_input ]]; then
        DEAL_II_DIR=$deal_ii_dir_input
    fi
    if [[ $number_jobs_input =~ ^[0-9]+$ ]]; then
        NUMBER_JOBS=$number_jobs_input
    fi
    if [[ $run_test_input =~ ^[Yy]$ ]]; then
        RUN_TEST=true
    fi
}

function print_configuration {
    echo ""
    echo "Configuration to be installed"
    echo "-----------------------------------"
    if [ $INSTALL_ASSIMP == true ]; then
        echo "Install assimp-v$ASSIMP_VERSION in \"$ASSIMP_DIR\""
    else
        echo "Use assimp installation in \"$ASSIMP_DIR\""
    fi
    if [ $INSTALL_PETSC == true ]; then
        echo "Install PETSc-v$PETSC_VERSION in \"$PETSC_DIR\" using PETSC_ARCH=$PETSC_ARCH and PETSC_FLAGS=\"$PETSC_FLAGS\""
    else
        echo "Use PETSc installation in \"$PETSC_DIR\""
    fi
    if [ $INSTALL_P4EST == true ]; then
        echo "Install p4est-v$P4EST_VERSION in \"$P4EST_DIR\""
    else
        echo "Use p4est installation in \"$P4EST_DIR\""
    fi
    if [ $INSTALL_DEAL_II == true ]; then
        echo "Install deal.II-v$DEAL_II_VERSION in \"$DEAL_II_DIR\" using DEAL_II_FLAGS=\"$DEAL_II_FLAGS\""
    else
        echo "Use deal.II installation in \"$DEAL_II_DIR\""
    fi
    echo "Compile using $NUMBER_JOBS job(s)"
    if [ $RUN_TEST == true ]; then
        echo "Run test suite after installation"
    else
        echo "Do not run test suite after installation"
    fi
    echo ""

    if [ $INTERACTIVE == true ]; then
        read -rp "Is this configuration correct? [y/N]: " confirm_input
        if [[ ! $confirm_input =~ ^[Yy]$ ]]; then
            echo "Installation aborted."
            exit 1
        fi
    fi
}

function install_prerequisites {
    echo "Are the following prerequisites installed on your machine?"
    echo ""
    echo "  $PREREQUISITES"
    echo ""
    read -rp "[y/N]: " install_prerequisite_input

    if [[ ! $install_prerequisite_input =~ ^[Nn]$ ]]; then
        return
    fi

    read -rp "Which package manager do you want to use to install the prerequisites? [brew/apt/apt-get/...]: " package_installer_input
    echo ""
    if [[ $package_installer_input ]]; then
        PACKAGE_INSTALLER=$package_installer_input
    fi

    if [ "$PACKAGE_INSTALLER" == "brew" ]; then
        echo "To install the prerequisites with brew, run the following command:"
        echo ""
        echo "  brew install $PREREQUISITES_BREW"
        echo ""
        echo "We recommend to pin the installed packages:"
        echo ""
        echo "  brew pin $PREREQUISITES_BREW"
        echo ""
    elif [ "$PACKAGE_INSTALLER" == "apt" ]; then
        echo "To install the prerequisites with apt, run the following command:"
        echo ""
        echo "  sudo apt install $PREREQUISITES_APT"
        echo ""
    elif [ "$PACKAGE_INSTALLER" == "apt-get" ]; then
        echo "To install the prerequisites with apt-get, run the following command:"
        echo ""
        echo "  sudo apt-get install $PREREQUISITES_APT"
        echo ""
    else
        echo "Please install the prerequisites manually."
        echo "Hint: If you are using a Debian-based distribution like Ubuntu, the following command can be used:"
        echo ""
        echo "  sudo apt install $PREREQUISITES_APT"
    fi
    echo "Afterwards you can continue with the installation with:"
    echo ""
    echo "  ./install_dealii.sh --skip-prerequisites"
    exit 2
}

function install_assimp {
    echo ""
    echo "Installing assimp-v$ASSIMP_VERSION in \"$ASSIMP_DIR\""
    echo "-----------------------------------"
    echo ""

    cd "$WORKPWD"
    tarname="v$ASSIMP_VERSION"
    dirname="assimp-$ASSIMP_VERSION"
    curl -LO "https://github.com/assimp/assimp/archive/refs/tags/$tarname.tar.gz"
    tar -xzf "$tarname.tar.gz"
    rm "$tarname.tar.gz"
    cd "$dirname"

    mkdir -p "$ASSIMP_DIR"
    rm -rf "$ASSIMP_DIR"
    cmake -S . -B build -DCMAKE_INSTALL_PREFIX="$ASSIMP_DIR"
    cd build
    make -j"$NUMBER_JOBS" install

    cd "$WORKPWD"
    rm -rf "$dirname"
}

function install_petsc {
    echo ""
    echo "Installing PETSc-v$PETSC_VERSION in \"$PETSC_DIR\" using PETSC_ARCH=$PETSC_ARCH and PETSC_FLAGS=\"$PETSC_FLAGS\""
    echo "-----------------------------------"
    echo ""

    cd "$WORKPWD"
    dirname="petsc-v$PETSC_VERSION"
    curl -LO "https://gitlab.com/petsc/petsc/-/archive/v$PETSC_VERSION/$dirname.tar.gz"
    tar -xzf "$dirname.tar.gz"
    rm "$dirname.tar.gz"
    mkdir -p "$PETSC_DIR"
    rm -rf "$PETSC_DIR"
    cd "$dirname"

    # shellcheck disable=SC2086
    ./configure --prefix="$PETSC_DIR" --PETSC_DIR="$(pwd)" --PETSC_ARCH="$PETSC_ARCH" $PETSC_FLAGS
    make -j"$NUMBER_JOBS" PETSC_DIR="$(pwd)" PETSC_ARCH="$PETSC_ARCH" all
    make -j"$NUMBER_JOBS" PETSC_DIR="$(pwd)" PETSC_ARCH="$PETSC_ARCH" install
    
    if [[ $RUN_TEST == true ]]; then
        make PETSC_DIR="$(pwd)" PETSC_ARCH="$PETSC_ARCH" check
        make PETSC_DIR="$PETSC_DIR" PETSC_ARCH="" check
    fi

    cd "$WORKPWD"
    rm -rf "$dirname"
}

function install_p4est {
    echo ""
    echo "Installing p4est-v$P4EST_VERSION in \"$P4EST_DIR\""
    echo "-----------------------------------"
    echo ""

    cd "$WORKPWD"
    dirname="p4est-$P4EST_VERSION"
    curl -LO "https://p4est.github.io/release/$dirname.tar.gz"
    curl -LO "https://dealii.org/current/external-libs/p4est-setup.sh"
    mkdir -p "$P4EST_DIR"
    rm -rf "$P4EST_DIR"
    bash p4est-setup.sh "$dirname.tar.gz" "$P4EST_DIR"

    rm -rf "$dirname"
    rm "$dirname.tar.gz"
    rm -rf p4est-build
    rm p4est-setup.sh

    cd "$WORKPWD"
}

function install_deal_ii {
    echo ""
    echo "Installing deal.II-v$DEAL_II_VERSION in \"$DEAL_II_DIR\" using DEAL_II_FLAGS=\"$DEAL_II_FLAGS\""
    echo "-----------------------------------"
    echo ""

    # Download and extract deal.II source code
    cd "$WORKPWD"
    dirname="dealii-$DEAL_II_VERSION"
    curl -LO "https://github.com/dealii/dealii/releases/download/v$DEAL_II_VERSION/$dirname.tar.gz"
    tar xf "$dirname.tar.gz"
    rm "$dirname.tar.gz"
    cd "$dirname"

    # Set CMake flags for deal.II configuration
    tmp_deal_ii_flags=(
        "-DCMAKE_INSTALL_PREFIX=$DEAL_II_DIR"
        "-DDEAL_II_CXX_FLAGS=-std=c++20"
        # "-DDEAL_II_WITH_TBB=ON"
        "-DDEAL_II_WITH_MPI=ON"
        "-DDEAL_II_WITH_HDF5=ON"
        "-DDEAL_II_WITH_ASSIMP=ON"
        "-DASSIMP_DIR=$ASSIMP_DIR"
        "-DDEAL_II_WITH_PETSC=ON"
        "-DPETSC_DIR=$PETSC_DIR"
        "-DPETSC_ARCH=$PETSC_ARCH"
        "-DDEAL_II_WITH_P4EST=ON"
        "-DP4EST_DIR=$P4EST_DIR"
    )

    # Print CMake flags for debugging purposes
    # shellcheck disable=SC2086
    echo Configure deal.II using "${tmp_deal_ii_flags[@]}" $DEAL_II_FLAGS

    # Configure and build deal.II
    mkdir -p "$DEAL_II_DIR"
    rm -rf "$DEAL_II_DIR"
    # shellcheck disable=SC2086
    cmake -S . -B build "${tmp_deal_ii_flags[@]}" $DEAL_II_FLAGS
    cd build
    make -j"$NUMBER_JOBS" install

    # Run test suite if requested
    if [ "$RUN_TEST" = true ]; then
        make test
    fi

    # Clean up
    cd "$WORKPWD"
    rm -rf "$dirname"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
    -n | --non-interactive)
        INTERACTIVE=false
        shift
        ;;
    -skp | --skip-prerequisites)
        SKIP_PREREQUISITE=true
        shift
        ;;
    -asp | --assimp)
        INSTALL_ASSIMP=true
        shift
        ;;
    -aspd | --assimp-dir)
        ASSIMP_DIR="$2"
        shift
        shift
        ;;
    -aspv | --assimp-version)
        ASSIMP_VERSION="$2"
        shift
        shift
        ;;
    -psc | --petsc)
        INSTALL_PETSC=true
        shift
        ;;
    -pscd | --petsc-dir)
        PETSC_DIR="$2"
        shift
        shift
        ;;
    -pscv | --petsc-version)
        PETSC_VERSION="$2"
        shift
        shift
        ;;
    -pscf | --petsc-flags)
        PETSC_FLAGS="$2"
        shift
        shift
        ;;
    -p4e | --p4est)
        INSTALL_P4EST=true
        shift
        ;;
    -p4ed | --p4est-dir)
        P4EST_DIR="$2"
        shift
        shift
        ;;
    -p4ev | --p4est-version)
        P4EST_VERSION="$2"
        shift
        shift
        ;;
    -d | --dealii)
        INSTALL_DEAL_II=true
        shift
        ;;
    -dd | --dealii-dir)
        DEAL_II_DIR="$2"
        shift
        shift
        ;;
    -dv | --dealii-version)
        DEAL_II_VERSION="$2"
        shift
        shift
        ;;
    -df | --dealii-flags)
        DEAL_II_FLAGS="$2"
        shift
        shift
        ;;
    -j | --jobs)
        NUMBER_JOBS="$2"
        shift
        shift
        ;;
    -t | --test)
        RUN_TEST=true
        shift
        ;;
    *)
        echo "Unknown option: $1"
        exit 1
        ;;
    esac
done

# Install prerequisites
if [ $SKIP_PREREQUISITE == false ] && [ $INTERACTIVE == true ]; then
    install_prerequisites
fi

# Set configuration interactively
if [ $INTERACTIVE == true ]; then
    set_configuration
fi

# Print configuration
print_configuration

echo "Starting installation. This may take a while!"

# Install assimp
if [ $INSTALL_ASSIMP == true ]; then
    install_assimp
fi

# Install PETSc
if [ $INSTALL_PETSC == true ]; then
    install_petsc
fi

# Install p4est
if [ $INSTALL_P4EST == true ]; then
    install_p4est
fi

# Install deal.II
if [ $INSTALL_DEAL_II == true ]; then
    install_deal_ii
fi

echo ""
echo "-----------------------------------"
echo "Installation complete."

exit 0
