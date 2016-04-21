fMain = open('likelihood.tex', 'w')

# Insert start of section
fMain.write("\input{likeStart.tex}\n\n\subsection{List of likelihoods}\n")

# Get list of likelihoods
from os import listdir
from os.path import isfile, join
likDir = "../../../../r-inla.org/doc/likelihood/"
likList = [f for f in listdir(likDir) if f.endswith('.tex')]

# Extract likelihood names
likNames = likList
for idx in range(0, len(likList)):
	parts = likList[idx].split('.')
	likNames[idx] = parts[0]
	
	# Add code for each likelihood
	fMain.write('\input{%s}\n' % join(likDir, likList[idx]))
	

# Close
fMain.close()
