export PATH=$HOME/p/inla/build/windows/bin:$PATH
if [ "${PS1:-}" != "" ]; then
    if [ "${WINDOW:-}" = "" ]; then
        PS1="MinGW[\W]\$ "
    fi
fi





