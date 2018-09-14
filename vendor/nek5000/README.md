# Nek5000 (cmake README)

## Overview

This branch of Nek5000 allows CMake to build the following components:
* Nek5000 library and driver
* gslib
* genmap
* genbox

The following components are **not** yet supported through CMake, but can be built through the existing build system:
* CVODE
* All Nek5000 tools other than genmap and genbox.  

Using CMake required the addition of CMakeList.txt files in the project subdirectories.  No changes to the existing 
source code were required.  The CMake build system exists alongside the existing build systems.  

## Usage

### Building Nek5000

The typical Nek5000 workflow involves configuring and building Nek5000 in a directory for a specific problem, which 
contains the SIZE, .usr, and mesh files needed for the problem.  This is still the general workflow with CMake.  For
example, to build eddy_uv:
``` Console
$ export CC=mpicc
$ export FC=mpif90
$ cd examples/eddy_uv
$ cmake -DCASENAME=eddy_uv ../../
$ make
```
The above commands illustrate new requirements of the CMake build system:
*. CMake uses the *environment variables* `CC` and `FC` to detect the compilers.  In this example, they are set using 
   the `export` commands.  
*. The casename is specified by the *CMake variable* `CASENAME`.  These are passed to CMake using the 
   `-D<variable>=<value>` syntax.
   
The resulting executables and libraries are:
* `core/nek5000`
* `core/libnek.a`
* `3rd_party/gslib/gslib-1.0.1/libgs.a`

The executable `core/nek5000` can be run as normal.  For example:
``` Console
$ mpirun -np 4 core/nek5000
```
Note that the `SESSION.NAME` file is created during the configure phase and need not be created at runtime.  


### Building tools

After configuring (i.e., calling `cmake -DCASENAME=eddy_uv ../../`), you may also compile genmap or genmap by invoking 
make:
``` Console
$ make genmap
$ make genbox
```
This produces `tools/genmap` and `tools/genbox`.  The usage of the tools are the same.  

### Building gslib standalone

gslib may also be build standalone, using the CMakeLists.txt in `3rd_party/gslib/gslib-1.0.1`.


# Nek5000 (master README)

| **`Short Tests`** | **`Examples`** |
|-----------------|---------------------|
| [![Build](https://travis-ci.org/Nek5000/Nek5000.svg?branch=master)](https://travis-ci.org/Nek5000/Nek5000) | [![Build Status](https://jenkins-ci.cels.anl.gov/buildStatus/icon?job=Nek5000)](https://jenkins-ci.cels.anl.gov/job/Nek5000/) |

Nek5000 is a fast and scalable open source CFD solver.

## Release Notes
Make sure to read the [release notes](https://github.com/Nek5000/Nek5000/blob/master/RELEASE.md) before using the code.


## Documentation

Visit our  [User's Guide](http://Nek5000.github.io/NekDoc/).

## Troubleshooting

If you run into problems compiling, installing, or running Nek5000, first check the User's Guide. If you are not able to find a solution to your problem there, please send a message to the User's Group [mailing list](https://lists.mcs.anl.gov/mailman/listinfo/nek5000-users).

## Reporting Bugs
Nek5000 is hosted on GitHub and all bugs are reported and tracked through the [Issues](https://github.com/Nek5000/Nek5000/issues) feature on GitHub. However, GitHub Issues should not be used for common troubleshooting purposes. If you are having trouble installing the code or getting your model to run properly, you should first send a message to the User's Group mailing list. If it turns out your issue really is a bug in the code, an issue will then be created on GitHub. If you want to request that a feature be added to the code, you may create an Issue on GitHub.

## Contributing

Our project is hosted on [GitHub](https://github.com/Nek5000). If you are planning a large contribution, we encourage you to discuss the concept here on GitHub and interact with us frequently to ensure that your effort is well-directed.

### How we do it
- Anything in master is always deployable
- Upcoming feature release get their own tags or branch that are branched out of master
- All development happens on the master branch.
- To work on something new, create a short lived local branch off of master
- When you need feedback or help, or your change is ready for merging, open a pull request.

### One-time Setup
1. Fork our [GitHub project](https://github.com/Nek5000/Nek5000)
2. Download fork with `git clone -o myfork https://github.com/<username>/Nek5000.git ~/Nek5000`
3. Add our repo `cd ~/Nek5000; git remote add origin https://github.com/Nek5000/Nek5000.git`
4. Download our repo `git fetch origin`
5. Set upstream for local master branch `git branch --set-upstream master remotes/origin/master`
6. Run `~/Nek5000/bin/git-hub setup —u <your username on GitHub> --global`
7. Add this to your [hub] section in `~/.gitconfig`:

```
[hub]
        ...
        upstream = Nek5000/Nek5000
        forkremote = myfork
```

### Typical Workflow
1. Create a feature branch hosting your change with `nekgit_co <descriptive name>`. Using a dedicated branch for every feature helps you to move between different developments while some are work in progress or under review.
2. Implement your code changes. To reset your branch and discard any changes run `git reset --hard origin/master`. To revert a set of files run `git checkout file1 file2 ...`
3. Commit your changes to your local repo using e.g. `git commit file1 file2 ...`. Do this frequently to save your work.
4. Periodically, changes made in our master should be pulled back into your local branch by `git pull -r`. This ensures that we do not end up in integration hell that will happen when many feature branches need to be combined at once.
5. If there are no merge conflicts, go to the next step. In case of conflicts edit the unmerged files in question. Merge conflicts are indicated  by the conflict marker `<<<<<<<` in your file.
6. Assuming you are happy run `nekgit_push`. This will create a pull request on GitHub. You can check with `git diff origin/master` what your push will do. When your pull request was merged, run `git pull` on your local master branch to see your change. You can delete the branch created in step (1) with `nekgit_rm <my branch name>`.
