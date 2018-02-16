# ParaFEM-lite

### What is this?

This is a streamlined version of the [ParaFEM](http://parafem.org.uk) codebase currently located on
sourceforge.

The original custom build system has been replaced with CMake, and the old bundled linear algebra
libraries removed in favour of modern libraries optimised for current hardware.


### Building

Building is done using an "out-of-tree" CMake build.  This just means you don't build right in the 
source folder.  Just make a new build folder and run cmake there:

```
git clone git@github.com:ptooley/ParaFEM-lite.git
cd ParaFEM-lite
mkdir build
cd build
cmake ..
make
```

### Status

The new build system is incomplete so far only building what I immediately need.  I intend to
come back and finish the job at some point though.

### Support

I did this for my own sanity because I couldn't get the original build scripts to work properly for
me. I never really set out to support a repackaged version of ParaFEM. That said if this is useful
to you please drop me a message, and if you have a problem please open an issue.  I can't promise I
will get back to you but I will try.
