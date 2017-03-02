for file in figs/*.pdf
do
  pdftops -eps "${file}"
done
