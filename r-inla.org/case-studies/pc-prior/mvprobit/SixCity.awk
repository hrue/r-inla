BEGIN {
    print "smoking", "age7", "age8", "age9", "age10";
}
NF == 0 {
    next;
}
/^smoking/{
    smoking = $2;
    next;
}
{
    n = $1;
    age7 = $2;
    age8 = $3;
    age9 = $4;
    age10 = $5;
    for(i = 1;  i<= n;  i=i+1) {
	print smoking, age7, age8, age9, age10;
    }
}


	
    
	
    
