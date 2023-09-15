#!/bin/bash
set -e

# Set default values
INTERACTIVE=false
INSTALL_PREREQUISITE=false
if [[ $(uname) == "Darwin" ]]; then
    PACKAGE_INSTALLER="brew"
else
    PACKAGE_INSTALLER="apt"
fi

INSTALL_PETSC=false
PETSC_DIR=${PETSC_DIR:-"$HOME/.local/lib/petsc"}
PETSC_VERSION="3.19.5"
if [[ $(uname) == "Darwin" ]]; then
    PETSC_ARCH=${PETSC_ARCH:-"arch-darwin-c-debug"} #TODO: Use debug or opt: arch-darwin-c-opt?
else
    PETSC_ARCH=${PETSC_ARCH:-"arch-linux-c-debug"}
fi

INSTALL_SLEPC=false
SLEPC_DIR=${SLEPC_DIR:-"$HOME/.local/lib/slepc"}
SLEPC_VERSION="3.19.2"

INSTALL_P4EST=false
P4EST_DIR=${P4EST_DIR:-"$HOME/.local/lib/p4est"}
P4EST_VERSION="2.8.5"

INSTALL_DEAL_II=false
DEAL_II_DIR=${DEAL_II_DIR:-"$HOME/.local/lib/dealii"}
DEAL_II_VERSION="9.5.1"

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
    read -p "Do you want to install prerequisites? [y/N]: " install_prerequisite_input
    read -p "Do you want to install PETSc? [y/N]: " install_petsc_input
    read -p "Enter the installation directory for PETSc (default: $PETSC_DIR): " petsc_dir_input
    read -p "Do you want to install SLEPc? [y/N]: " install_slepc_input
    read -p "Enter the installation directory for SLEPc (default: $SLEPC_DIR): " slepc_dir_input
    read -p "Do you want to install p4est? [y/N]: " install_p4est_input
    read -p "Enter the installation directory for p4est (default: $P4EST_DIR): " p4est_dir_input
    read -p "Do you want to install deal.II? [y/N]: " install_deal_ii_input
    read -p "Enter the installation directory for deal.II (default: $DEAL_II_DIR): " deal_ii_dir_input
    read -p "Enter the number of processors to use for compilation (default: $NUMBER_JOBS): " number_jobs_input
    read -p "Do you want to run the test suite after installation? [y/N]: " run_test_input

    # Set the configuration based on user input
    if [[ $install_prerequisite_input =~ ^[Yy]$ ]]; then
        INSTALL_PREREQUISITE=true
    fi
    if [[ $install_petsc_input =~ ^[Yy]$ ]]; then
        INSTALL_PETSC=true
    fi
    if [[ $petsc_dir_input ]]; then
        PETSC_DIR=$petsc_dir_input
    fi
    if [[ $install_slepc_input =~ ^[Yy]$ ]]; then
        INSTALL_SLEPC=true
    fi
    if [[ $slepc_dir_input ]]; then
        SLEPC_DIR=$slepc_dir_input
    fi
    if [[ $install_p4est_input =~ ^[Yy]$ ]]; then
        INSTALL_P4EST=true
    fi
    if [[ $p4est_dir_input ]]; then
        P4EST_DIR=$p4est_dir_input
    fi
    if [[ $install_deal_ii_input =~ ^[Yy]$ ]]; then
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
    if [ $INSTALL_PREREQUISITE == true ] ; then
        echo "Install prerequisites using $PACKAGE_INSTALLER"
    else
        echo "Not installing prerequisites"
    fi
    if [ $INSTALL_PETSC == true ] ; then
        echo "Install PETSc-v$PETSC_VERSION in \"$PETSC_DIR\""
    else
        echo "Use PETSc installation in $PETSC_DIR"
    fi
    if [ $INSTALL_SLEPC == true ] ; then
        echo "Install SLEPc-v$SLEPC_VERSION in \"$SLEPC_DIR\""
    else
        echo "Use SLEPc installation in $SLEPC_DIR"
    fi
    if [ $INSTALL_P4EST == true ] ; then
        echo "Install p4est-v$P4EST_VERSION in \"$P4EST_DIR\""
    else
        echo "Use p4est installation in $P4EST_DIR"
    fi
    if [ $INSTALL_DEAL_II == true ] ; then
        echo "Install deal.II-v$DEAL_II_VERSION in \"$DEAL_II_DIR\""
    else
        echo "Use deal.II installation in $DEAL_II_DIR"
    fi
    echo "Compile on $NUMBER_JOBS processor(s)"
    echo "RUN_TEST=$RUN_TEST"
    echo ""
}

#TODO: Validate configuration (e.g. check if directories exist)

function install_prerequisites {
    echo "Installing prerequisites using $PACKAGE_INSTALLER"
    if [ $PACKAGE_INSTALLER == "brew" ] ; then
        brew install gcc make cmake open-mpi hdf5-mpi boost lapack tbb
    elif [ $PACKAGE_INSTALLER == "apt" ] ; then
        sudo apt install gcc make cmake openmpi-bin libhdf5-openmpi-dev libboost-all-dev liblapack-dev libtbb2 libtbb2-dev
    else
        echo "$PACKAGE_INSTALLER is not a supported package installer. Use one of brew, apt."
        exit 1
    fi
}

function install_petsc {
    echo ""
    echo "Installing PETSc-v$PETSC_VERSION in \"$PETSC_DIR\" using PETSC_ARCH=$PETSC_ARCH"
    echo "-----------------------------------"
    echo ""

    cd "$WORKPWD"
    dirname="petsc-v$PETSC_VERSION"
    curl -LO "https://gitlab.com/petsc/petsc/-/archive/v$PETSC_VERSION/$dirname.tar.gz"
    tar -xzf "$dirname.tar.gz"
    rm "$dirname.tar.gz"
    rm -rf "$PETSC_DIR"
    mv "$dirname" "$PETSC_DIR"
    cd "$PETSC_DIR"

    ./configure --download-hypre --download-scalapack --download-mumps
    make -j"$NUMBER_JOBS" PETSC_DIR="$PETSC_DIR" PETSC_ARCH="$PETSC_ARCH" all

    if [[ $RUN_TEST == true ]]; then
        make PETSC_DIR="$PETSC_DIR" PETSC_ARCH="$PETSC_ARCH" check
    fi

    cd "$WORKPWD"
}

function install_slepc {
    echo ""
    echo "Installing SLEPc-v$SLEPC_VERSION in \"$SLEPC_DIR\""
    echo "-----------------------------------"
    echo ""

    cd "$WORKPWD"
    dirname="slepc-$SLEPC_VERSION"
    curl -LO "https://slepc.upv.es/download/distrib/$dirname.tar.gz"
    tar -xzf "$dirname.tar.gz"
    rm "$dirname.tar.gz"
    rm -rf "$SLEPC_DIR"
    mv "$dirname" "$SLEPC_DIR"
    cd "$SLEPC_DIR"

    ./configure
    make -j"$NUMBER_JOBS" SLEPC_DIR="$SLEPC_DIR" PETSC_DIR="$PETSC_DIR" PETSC_ARCH="$PETSC_ARCH" all

    if [[ $RUN_TEST == true ]]; then
        make SLEPC_DIR="$SLEPC_DIR" PETSC_DIR="$PETSC_DIR" PETSC_ARCH="$PETSC_ARCH" check
    fi

    cd "$WORKPWD"
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
    echo "Installing deal.II-v$DEAL_II_VERSION in \"$DEAL_II_DIR\""
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
    deal_ii_flags=(
        "-DCMAKE_INSTALL_PREFIX=$DEAL_II_DIR"
        "-DDEAL_II_CXX_FLAGS=-std=c++17"
        "-DDEAL_II_WITH_MPI=ON"
        "-DDEAL_II_WITH_HDF5=ON"
        "-DDEAL_II_WITH_PETSC=ON"
        "-DPETSC_DIR=$PETSC_DIR"
        "-DPETSC_ARCH=$PETSC_ARCH"
        "-DDEAL_II_WITH_SLEPC=ON"
        "-DSLEPC_DIR=$SLEPC_DIR"
        "-DDEAL_II_WITH_P4EST=ON"
        "-DP4EST_DIR=$P4EST_DIR"
    )

    # Print CMake flags for debugging purposes
    echo "${deal_ii_flags[@]}"

    # Configure and build deal.II
    rm -rf "$DEAL_II_DIR"
    cmake -S . -B build  "${deal_ii_flags[@]}"
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
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
        -i|--interactive)
        INTERACTIVE=true
        shift
        ;;
        -p|--prerequisites)
        INSTALL_PREREQUISITE=true
        shift
        ;;
        -pm|--package-manager)
        PACKAGE_INSTALLER="$2"
        shift
        shift
        ;;
        -psc|--petsc)
        INSTALL_PETSC=true
        shift
        ;;
        -pscd|--petsc-dir)
        PETSC_DIR="$2"
        shift
        shift
        ;;
        -pscv|--petsc-version)
        PETSC_VERSION="$2"
        shift
        shift
        ;;
        -slc|--slepc)
        INSTALL_SLEPC=true
        shift
        ;;
        -slcd|--slepc-dir)
        SLEPC_DIR="$2"
        shift
        shift
        ;;
        -slcv|--slepc-version)
        SLEPC_VERSION="$2"
        shift
        shift
        ;;
        -p4e|--p4est)
        INSTALL_P4EST=true
        shift
        ;;
        -p4ed|--p4est-dir)
        P4EST_DIR="$2"
        shift
        shift
        ;;
        -p4ev|--p4est-version)
        P4EST_VERSION="$2"
        shift
        shift
        ;;
        -d|--dealii)
        INSTALL_DEAL_II=true
        shift
        ;;
        -dd|--dealii-dir)
        DEAL_II_DIR="$2"
        shift
        shift
        ;;
        -dv|--dealii-version)
        DEAL_II_VERSION="$2"
        shift
        shift
        ;;
        -np|--number-of-processors)
        NUMBER_JOBS="$2"
        shift
        shift
        ;;
        -t|--test)
        RUN_TEST=true
        shift
        ;;
        *)
        echo "Unknown option: $1"
        exit 1
        ;;
    esac
done

# Set configuration interactively
if [ $INTERACTIVE == true ] ; then
    set_configuration
fi

# Print configuration
print_configuration

# Ask for confirmation
read -p "Do you want to continue with the installation? [y/N]: " confirm_input
# Abort the installation if the user does not confirm
if [[ ! $confirm_input =~ ^[Yy]$ ]]; then
    echo "Installation aborted."
    exit 1
fi

# Install prerequisites
if [ $INSTALL_PREREQUISITE == true ] ; then
    install_prerequisites
fi

# Install PETSc
if [ $INSTALL_PETSC == true ] ; then
    install_petsc
fi

# Install SLEPc
if [ $INSTALL_SLEPC == true ] ; then
    install_slepc
fi

# Install p4est
if [ $INSTALL_P4EST == true ] ; then
    install_p4est
fi

# Install deal.II
if [ $INSTALL_DEAL_II == true ] ; then
    install_deal_ii
fi

echo ""
echo "-----------------------------------"
echo "Installation complete."

exit 0