fp  = fopen('design-f.dat','w')

for i=2:11
    
    D = ccdesign(i,'type', 'f')
    D = unique(D,'rows');
    
    fprintf(fp,'%d %d\n', size(D'));
    m=size(D);
    for jj=1:m(1)
	for kk=1:m(2)
	    fprintf(fp,' %.10f', D(jj,kk));
	end
	fprintf(fp,'\n');
    end
    
end
fclose(fp)

