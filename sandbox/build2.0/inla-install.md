** Outline of the new R-INLA install procedure **

This is an outline of the proposed R-INLA install procedure, designed
to

- separate the R code from the precompiled binaries (which potentially
  can make it into CRAN in the future)
- make it easy for the normal user to install the inla-program
- allow for sufficient flexibility for the expert-user/developer

** Dependency of ```libR``` and ```R-version``` **

The new build system does not depend on the version of ```R``` and/or
```libR```. This because the inla-program now read the ```libR```
corresponding to the version of ```R``` it being called from as the
env-variable ```R_HOME``` is set automatically, which gives the path
to the ```libR```. If the inla-program is lauched outside of ```R```,
like the terminal, then we need to set ```R_HOME``` manually

```
env R_HOME=$(R RHOME) inla-program -v -t8:1:2 Model.ini
```

This is required iff any of the ```libR``` spesific
features are used: ```rgeneric``` and ```rprior```.

** First time install and later updates (normal user) **

To install R-INLA for the first time, we first install the package
** and ** then the inla-program. To install the R-INLA package we do

```
pak("r-inla/install")
```

To install the inla-program, we do

```
INLA::inla.stiles.install()
```

which will install the inla-program matching the R-INLA version. The
inla-program is stored 

- within the R-INLA package itself (if there is write access),

- otherwise it is stored at a predefined location at the user's home
  directory, f.ex ```~/.cache/R/INLA/Version_xx.xx.xx/...``` for Linux
  and Mac

- if none of the two locations above have with write access, then exit
  with an error saying developer-install is required (described later)

This approach avoids storing ```inla.call``` in the global options,
which we could avoid since ```inla.call``` is not really an
``option''. Note there are no other changes required (like adding
setup-code to the user's ```.Rprofile```).

For later updates/upgrades, then we can follow the style used earlier,
and do one of

```
inla.upgrade()  
inla.upgrade(testing=TRUE)
```
This will will either upgrade the stable version or the most recent (testing)
one. The inla-program will be also be upgraded when calling
```inla.upgrade()```.

We have to check if a package install with ```pak()``` will be updated
using one of the standard upgrade functions already in ```R``` or in
one of common ```pak```-like packages. If this is the case, then the
inla-program will need to be installed again using
```inla.stiles.install()```.

In the case where the inla-program is installed within the R-INLA
package itself, the old inla-program will automatically be removed when the
R-INLA package is updated.

In the case where the inla-program is installed at a the user's home
directory, then the new inla-program will be added and but the old
inla-program is kept. This because the inla-program is stored with
path ```.../Version_xx.yy.zz/...'' hence the path is unique.
We could add a function to list installed inla-programs, delete
unused/older ones, etc, to avoid things to pile up.

** First time install and later updates (developer/expert) **

This pathway is target towards expert-users/developers, and will
mainly be used for testing and development purposes.

The normal procedure, is to install ```R-INLA``` and the ```inla-program```
as a normal-user. 

This pathway add extentions to

- install and run with another version of the inla-program
- use a user-compiled version of the inla-program. 

We need to add a function to list the available inla-programs
available for download. One of those can then be given as argument to
```inla.stiles.install(version="...")```. This inla-program is installed
along-side of the original one, again with path
```.../Version_xx.yy.zz/...``` so the path is unique.

To allow for multiple running inla-programs with different versions,
we do not allow to store (as we did before) the path in the global 
```INLA options```.  The path to the inla-program needs to be 
given as the ```inla.call```
argument (that existed before as well)

``` 
inla(..., inla.call="/home/hrue/inla/devel/bin/inla.run")
inla(..., inla.call="Version_26.09.03") 
inla(..., inla.call="26.09.03") 
```

However, this approach can be annoying as one might need to add that argument
to existing code manually. We can bypass this, if we also check and read the environment
variable ```INLA_CALL```, like

``` 
Sys.setenv(INLA_CALL="/home/hrue/inla/devel/bin/inla.run")
Sys.setenv(INLA_CALL="Version_26.09.03")
Sys.setenv(INLA_CALL="26.09.03") 
inla(...) 
```

If different versions are running simultanously within the same
R-session, like in a parallel loop/regions, then the environment
variable ```INLA_CALL``` needs to be set within the parallel region
where it is thread private and all is ok.

If the ```inla.call```-argument or ```INLA_CALL``` environment
variable is used, then there should be no check if the versions of the
inla-program and R-INLA is the same. If they are not used, the
```inla()```-call will fail if a inla-program is not found with same
version.

** Note** Running older inla-programs with newer R-INLA versions, is not
guaranteed to work, as there could be additional entries in
```Model.ini``` file or featues used that are not there in older
inla-programs. Running newer inla-programs with older R-INLA versions,
will likely be more successful due to reasonable default values for
new options added, but there is no guarantee (although we will try to
make it backward compatible when this can be done easily).
