  for (i=0; i<(int)N[0] ; i++) {
/*     do{*/
        sobseqmprof(nf,x);
	X=Img[0] + (Img[mrows -1] - Img[0])*x[1];
/*      X=Img[mrows*0 + 0]+ (Img[mrows -1] -Img[mrows*0 + 0])/mrow*x[0]; */
	j=0;
/*
	do {

  	    printf ("%e \t %e \t %d \t %e \n", x[1], X, index, Img[j]);  

	    index= j ;
	    j++;
	    
	} while ((X<Img[j]));
	index=index;
	fX = Img[index] + (Img[mrows + index + 1] - Img[mrows + index])/ 
	        (Img[index + 1] - Img[index])* (X-Img[index]);
	rnd=x[1];
     } while ((rnd>fX));
*/
     XX[inc]=X;
     inc++;
	
  } 
