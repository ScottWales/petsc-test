#!/bin/bash
## \file    bootstrap.sh
#  \author  Scott Wales <scott.wales@unimelb.edu.au>
#
#  Copyright 2014 ARC Centre of Excellence for Climate Systems Science
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

# Script for Centos VM setup

# CERN development
wget -nv -O /etc/yum.repos.d/slc6-devtoolset.repo http://linuxsoft.cern.ch/cern/devtoolset/slc6-devtoolset.repo

# EPEL
rpm -i http://dl.fedoraproject.org/pub/epel/6/i386/epel-release-6-8.noarch.rpm

cat > /etc/pki/rpm-gpg/RPM-GPG-KEY-cern << EOF
-----BEGIN PGP PUBLIC KEY BLOCK-----
Version: GnuPG v1.2.6 (GNU/Linux)

mQGiBEK/0MURBACv5Rm/jRnrbyocW5t43hrjFxlw/DPLTWiA16apk3P2HQQ8F6cs
EY/gmNmUf4U8KB6ncxdye/ostSBFJmVYh0YEYUxBSYM6ZFui3teVRxxXqN921jU2
GbbWGqqlxbDqvBxDEG95pA9oSiFYalVfjxVv0hrcrAHQDW5DL2b8l48kGwCgnxs1
iO7Z/5KRalKSJqKx70TVIUkD/2YkkHjcwp4Nt1pPlKxLaFp41cnCEGMEZVsNIQuJ
1SgHyMHKBzMWkD7QHqAeW3Sa9CDAJKoVPHZK99puF8etyUpC/HfmOIF6jwGpfG5A
S7YbqHX6vitRlQt1b1aq5K83J8Y0+8WmjZmCQY6+y2KHOPP+zHWKe5TJDeqDnN0j
sZsKA/441IF4JJTPEhvRFsPJO5WKg1zGFbxRPKvgi7+YY6pJ0VFbOMcJVMkvSZ2w
4QRD+2ets/pRxNhITHfPToMV3lhC8m1Je5fzoSvSixgH/5o9mekWWSW7Uq7U0IWA
7OD7RraJRrGxy0Tz3G+exA7svv/zn9TW/BaHFlMHoyyDHOYZmIhhBB8RAgAhBQJC
v+/uFwyAEeb+6rc8Txi4s8pfgZAf4xOTel99AgcAAAoJEF4D/eUdHgNLGCgAmwdu
KegSOBXpDe061zF8NoN6+OFiAJ9nKo+uC6xBZ9Ey550SmhFCPPA2/rRTQ0VSTiBM
aW51eCBTdXBwb3J0IChSUE0gc2lnbmluZyBrZXkgZm9yIENFUk4gTGludXggU3Vw
cG9ydCkgPGxpbnV4LnN1cHBvcnRAY2Vybi5jaD6IWgQTEQIAGgUCQr/QxQULBwMC
AQMVAgMDFgIBAh4BAheAAAoJEF4D/eUdHgNL/HsAn1ntKwRoSA9L0r8UyF7Zqn3U
79m1AJ9Y4NsSE/dlFYLfmf0+baoq7b5asIicBBMBAgAGBQJCv9DjAAoJEPy9YCiW
u335GTwEALjUQ7+cHxi0sifstCLoyRYQSu7Eas0M1UD2ZxSQNBnYsx4rDZJk9TmK
7QCzR1yRw9aixzZsRlNbed5VPxSzn89PE5m7Sx1bpl89sPgZ4BY95AL2wExyDWRp
1ON2+ztYeYtT7ZCkmeM+PBzt6RHR/jo3361faBS+qOkmpiiRWf3XiEYEExECAAYF
AkK/0WAACgkQkB/jE5N6X33DFQCgkvy1ruogu5Ibs5CzGY/ALiSJhyAAn3ygn3p/
xrNQ8Dy5j4KfgJINoxT4iEYEEBECAAYFAkK/9CcACgkQDIloXtlLxZSiRACdG0kT
KlB4X4VBocUyxMReO9e5MvsAoIKWgcJYE8AGmRXjfIisCAzPtVX+iEYEExECAAYF
AkK/8oUACgkQtQgG0wyY/52z1ACgkkxNdhHKbEol3Kwka1tICWHMIwIAn3PWJQR0
C1MV1+CnT8UupHzxy6J7
=IUD3
-----END PGP PUBLIC KEY BLOCK-----
EOF
yum install --assumeyes devtoolset-2-toolchain openmpi-devel lapack-devel python-argparse

cat > /etc/profile.d/xxx-development.sh << EOF
# This needs to come after /etc/profile.d/modules.sh
source /opt/rh/devtoolset-2/enable
module load openmpi-x86_64

export LD_LIBRARY_PATH=/home/vagrant/petsc/lib:$LD_LIBARARY_PATH
EOF

ln -s /vagrant ~vagrant/petsc-test

wget -nv -O - http://walesnix.earthsci.unimelb.edu.au/petsc-3.4.4-debug.tar.gz | tar xz -C ~vagrant
wget -nv -O - http://walesnix.earthsci.unimelb.edu.au/pfunit-2.1.0.tar.gz | tar xz -C ~vagrant
