fMain = open('latent.tex', 'w')

# Insert start of section
fMain.write("\\input{latent/latentStart.tex}\cleardoublepage\n")

# Get list of likelihoods
from os import listdir
from os.path import isfile, join
likDir = "../../../../../r-inla.org/doc/latent/"
likList = [f for f in listdir(likDir) if f.endswith('.tex')]

# Extract likelihood names
for idx in range(0, len(likList)):
	# Get names of likelihoods
	currLik = likList[idx]
	parts = currLik.split('.')
	
#	# Add subsection for likelihood
#	fMain.write('\\subsubsection{\texttt{%s}}\n' % parts[0])

	# Read text from file
	fSource = open(join(likDir, likList[idx]), 'r')
	inPreamble = True
	inSection  = False
	for line in fSource:
		if 'end{document}' in line:
			if (inPreamble == True):
				fSource.close()
				break
			inPreamble = True
		if 'section*' in line:
			inSection = True
			inPreamble = True
		if (inSection == True):
			test = line.split('*')
			print(test)
			if(len(test) == 1):
				fMain.write(test[0])
			else:
				fMain.write(test[0]+test[1])
			if ('}' in line):
				inSection = False
				inPreamble = False
			continue
		if (inPreamble == False):
			if ('\\verbatiminput{' in line):
				if (line[0] == '{'):
					fMain.write('{')
					line = line[1:len(line)]
				test = line.split('{')
				fMain.write(test[0]+'{'+likDir+test[1])
			elif ('\\input{' in line):
				test = line.split('/')
				if ('hazard' in test[len(test)-2]):
					fMain.write('\\input{../../../../r-inla.org/doc/hyper/hazard/' + test[len(test)-1])
				elif ('likelihood' in test[len(test)-2]):
					fMain.write('\\input{../../../../r-inla.org/doc/hyper/likelihood/' + test[len(test)-1])
				elif ('latent' in test[len(test)-2]):
					fMain.write('\\input{../../../../r-inla.org/doc/hyper/latent/' + test[len(test)-1])
				else:
					fMain.write('\\input{../../../../r-inla.org/doc/hyper/prior/' + test[len(test)-1])
			else:
				fMain.write(line)
		
			
	fSource.close()

	# New page between each likelihood
	fMain.write('\\cleardoublepage\n')
	

# Close
fMain.close()
