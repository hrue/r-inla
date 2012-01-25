## keep track of updates needed to be done

all:; \
        ## update the doc extracted from the R-files
        -make -C r-inla.org/doc
        ## update the TAGS-file and the man-files that are generated
        ## automatically 
        -make -C rinla/R 
        ## update the TAGS-file
        -make -C gmrflib TAGS
        ## update the TAGS-file
        -make -C inlaprog src/TAGS