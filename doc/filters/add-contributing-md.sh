#!/bin/bash
# -----------------------------------------------------------------------------
#
# Copyright (C) 2023 by the Sapphire++ authors
#
# This file is part of Sapphire++.
#
# Sapphire++ is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Sapphire++ is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Sapphire++. If not, see <https://www.gnu.org/licenses/>.
#
# -----------------------------------------------------------------------------

# This filter avoids to double the CONTRIBUTING.md file. A part of it is in the
# root dir to make it visible on github.com. The rest is under doc/pages

# This script uses GNU sed. On MacOS, we use gsed, if available
if command -v gsed &>/dev/null; then
    SED=gsed
elif sed --version >/dev/null 2>&1; then
    SED=sed
else
    echo "!WARNING! This script requires GNU sed, but it's not installed. Please install it and try again."
    exit 1
fi

# The absolute path to the CONTRIBUTING.md file is
# given as a command-line argument (see Doxyfile.in)
CONTENT_CONTRIBUTING=$(cat $1)

# Since sed is sensible to a set of special characters, we have to escape these
# first.
CONTENT_ESCAPED=$($SED -e 's/[&\\/]/\\&/g; s/$/\\/' -e '$s/\\$//' <<<"${CONTENT_CONTRIBUTING}")

# Insert the content of CONTRIBUTING.md at the beginning of
# doc/pages/contributing.md and print to standard out
$SED "1i${CONTENT_ESCAPED}" "${2}"
