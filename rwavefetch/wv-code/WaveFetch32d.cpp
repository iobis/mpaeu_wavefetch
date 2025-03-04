/* Wave fetch models in C++ for R */

#include <Rcpp.h>
#include <cmath>
#include <RProgress.h>
// [[Rcpp::depends("progress")]]
using namespace Rcpp;

//' Compute the mean of a numeric vector.
//'
//' @param x NumericVector of values.
//' @return The arithmetic mean of x.
// [[Rcpp::export]]
double meanC(NumericVector x) {
  int n = x.size();
  double total = 0;
  
  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total / n;
}



//' Identify coastal cells in a land matrix.
//'
//' This function examines each cell and its 3×3 neighborhood in the provided land matrix.
//' A cell is considered coastal if it is land (value 1) but not completely surrounded by land.
//'
//' @param land IntegerMatrix representing land (1) and water (0).
//' @return An IntegerMatrix with coastal cells marked as 1.
// [[Rcpp::export]]
IntegerMatrix isitcoast(IntegerMatrix land) {
  int ncolx = land.ncol();
  int nrowx = land.nrow();
  //Rcout << "ncols" << ncolx << "nrows" << nrowx << std::endl;
  
  IntegerMatrix coastcell(nrowx,ncolx);
  
  for(int j = 0; j < ncolx; j++) {
    for(int i = 0; i < nrowx; i++) {
      int sumland = 0;
      int ncells = 0;
      for(int jshift = -1; jshift < 2; jshift++) {
        for(int ishift = -1; ishift < 2; ishift++) {
          if(j+jshift > -1 && j+jshift < ncolx && i+ishift> -1 && i+ishift<nrowx) {
            sumland += land(i+ishift,j+jshift);
            ncells += 1;
          }
        }
      }
      if(land(i,j)==1 && sumland<ncells) {
        coastcell(i,j) = 1;
      }
    }
    //Rcout << j << std::endl;
  }
  return coastcell;
}



//' Identify near-coast cells based on a land matrix and a search distance.
//'
//' For cells identified as coastal, this function marks surrounding water cells (value 0) 
//' within a distance dx as near-coast (value 2).
//'
//' @param land IntegerMatrix representing land (1) and water (0).
//' @param dx Integer specifying the search distance (in cells).
//' @return An IntegerMatrix with coastal cells marked as 1 and near-coast cells marked as 2.
// [[Rcpp::export]]
IntegerMatrix isitnearcoast(IntegerMatrix land, int dx) {
  int ncolx = land.ncol();
  int nrowx = land.nrow();
  //Rcout << "ncols" << ncolx << "nrows" << nrowx << std::endl;
  
  IntegerMatrix coastcell(nrowx,ncolx);
  coastcell = isitcoast(land);
  
  for(int j = 0; j < ncolx; j++) {
    for(int i = 0; i < nrowx; i++) {
      if (coastcell(i,j)==1) {
        for(int jshift = -dx; jshift < dx+1; jshift++) {
          for(int ishift = -dx; ishift < dx+1; ishift++) {
            if(j+jshift > -1 && j+jshift < ncolx && i+ishift> -1 && i+ishift<nrowx) {
              if(land(i+ishift,j+jshift)==0) {
                coastcell(i+ishift,j+jshift) = 2;
              }
            }
          }
        }
      }
    }
    //Rcout << j << std::endl;
  }
  return coastcell;
}



//' Compute the wave fetch metric for coastal cells (Variant 1).
//'
//' This function computes a wave fetch metric (wx32) for each cell in a land matrix using random
//' square-root distances and angles. It performs an inner search and alternative mapping when needed.
//'
//' @param land1 IntegerMatrix representing the primary land area.
//' @param land2 IntegerMatrix representing the secondary land area used for alternative mapping.
//' @param dx Integer specifying the search distance for near-coast detection.
//' @param dwx Integer controlling the random search distance magnitude.
//' @param cl Integer used to identify coastal cells.
//' @param nsearches Integer specifying the number of random search iterations.
//' @param land1llx Lower left x-coordinate for land1.
//' @param land1uly Upper left y-coordinate for land1.
//' @param land1cellszx Double representing the cell size in the x-direction for land1.
//' @param land1cellszy Double representing the cell size in the y-direction for land1.
//' @param land2llx Lower left x-coordinate for land2.
//' @param land2uly Upper left y-coordinate for land2.
//' @param land2cellszx Double representing the cell size in the x-direction for land2.
//' @param land2cellszy Double representing the cell size in the y-direction for land2.
//' @param verbose Integer print messages to track progress. Use 0 for no messages, 1 for progress bar or 2 for text messages.
//' @return A NumericMatrix containing the computed wave fetch metric for each cell.
// [[Rcpp::export]]
NumericMatrix coastwx32twomapllranf1(IntegerMatrix land1, IntegerMatrix land2, int dx, int dwx, int cl, int nsearches,
                                    double land1llx, double land1uly, double land1cellszx, double land1cellszy, 
                                    double land2llx, double land2uly, double land2cellszx, double land2cellszy,
                                    int verbose = 1) {
  int ncolx = land1.ncol();
  int nrowx = land1.nrow();
  int lncols =land2.ncol();
  int lnrows =land2.nrow();
  
  RProgress::RProgress pb("Processing [:bar] Done: :percent | ETA: :eta", ncolx);
  if (verbose == 2) Rcout << "ncols" << ncolx << "nrows" << nrowx << std::endl;

  double pi = 3.14159265358979323846;
//  int nsearches = 5000;
  double vdwx = sqrt(dwx)*1000; // changed from int
  
  IntegerMatrix coastcell(nrowx,ncolx);
  coastcell = isitnearcoast(land1,dx);
  
  NumericMatrix wx32(nrowx,ncolx);
  
  if (verbose == 1) pb.tick(0);
  for(int j = 0; j < ncolx; j++) {
    if (verbose == 1) pb.tick();
    for(int i = 0; i < nrowx; i++) {
      if (coastcell(i,j)==cl) {
        NumericVector minfetch(32);
        for(int k = 0; k < 32; k++) {
          minfetch(k)=dwx*land1cellszy;
        }
        // Random square root distance and angle
        for(int count = 1; count < nsearches+1;count+=1) {
          
          double vdist = ((double) rand() / RAND_MAX) * vdwx/1000;  // changed from int
          double ranangle = 2 * pi * ((double) rand() / RAND_MAX);
          // double vx = vdist*vdist*sin(vang);
          // double vy = vdist*vdist*cos(vang);
          // 
          // int iranangle = rand() % 10000;
          // int iranlogdist = rand() % 10000;
          
          double randist = vdist*vdist;
          
          //Rcout << ranangle << "  " << ranlogdist << "  " << randist << std::endl;
          
          int jshift = randist*sin(ranangle);
          int ishift = randist*cos(ranangle);
          if(j+jshift > -1 && j+jshift < ncolx && i+ishift> -1 && i+ishift<nrowx) {
            if(land1(i+ishift,j+jshift)==1) {
              double yishift=land1cellszy*ishift;
              double xjshift=land1cellszx*jshift;
              double disttoland =  sqrt((xjshift*xjshift)+(yishift*yishift)); 
              int angletoland = 32 * (pi + atan2(xjshift,yishift))/(2*pi); 
              if (angletoland>31) {angletoland -=32;}
              // Rcout << jshift << "  " << ishift << "  " << yishift << "   "<< xjshift << std::endl;
              //  Rcout << ranlogdist << "  " << ranangle << "  " << angletoland << "   "<< disttoland << std::endl;
              if (minfetch(angletoland)>disttoland ) {
                minfetch(angletoland)=disttoland;
              }
            }
          }
          else {
            // alt-map search code here
            int lj = (land1llx - land2llx + land1cellszx * (j + jshift)) / land2cellszx;
            int li = (land2uly - land1uly + land1cellszy * (i + ishift)) / land2cellszy;
            if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
              if(land2(li,lj)==1) {
                double yishift=land1cellszy*ishift;
                double xjshift=land1cellszx*jshift;
                double disttoland = sqrt((xjshift*xjshift)+(yishift*yishift)); 
                int angletoland = 32 * (pi + atan2(xjshift,yishift))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << jshift << "  " << ishift << "  " << yishift << "   "<< xjshift << std::endl;
                // Rcout << ranlogdist << "  " << ranangle << "  " << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
          }
        }
        
		        // Inner search         
        for(int jshift = -10; jshift < 11; jshift+=1) {
          for(int ishift = -10; ishift < 11; ishift+=1) {
            if(j+jshift > -1 && j+jshift < ncolx && i+ishift> -1 && i+ishift<nrowx) {
              if(land1(i+ishift,j+jshift)==1) {
                double yishiftj=land1cellszy*ishift;
                double xjshiftj=land1cellszx*jshift;
                double disttoland = sqrt((xjshiftj*xjshiftj)+(yishiftj*yishiftj)); 
                int angletoland = 32 * (pi + atan2(xjshiftj,yishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                if (angletoland<0) {angletoland +=32;}
                if (minfetch(angletoland)>disttoland ) {minfetch(angletoland)=disttoland;}
                // for (int ashift = -1; ashift<2; ashift +=1) {
                //   int angleref = angletoland + ashift;
                //   if (angleref>31) {angleref -=32;}
                //   if (angleref<0) {angleref +=32;}
                //   // Rcout << "angleref innerloop" << angleref << "   "<< disttoland << std::endl;
                //   if (minfetch(angleref)>disttoland ) {
                //     minfetch(angleref)=disttoland;
                //   }
                //}
              }
            } 
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellszx * (j + jshift)) / land2cellszx;
              int li = (land2uly - land1uly + land1cellszy * (i + ishift)) / land2cellszy;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double yishiftj=land1cellszy*ishift;
                  double xjshiftj=land1cellszx*jshift;
                  double disttoland = sqrt((xjshiftj*xjshiftj)+(yishiftj*yishiftj)); 
                  int angletoland = 32 * (pi + atan2(xjshiftj,yishiftj))/(2*pi); 
                  if (angletoland>31) {angletoland -=32;}
                  // Rcout << angletoland << "   "<< disttoland << std::endl;
                  if (minfetch(angletoland)>disttoland ) {
                    minfetch(angletoland)=disttoland;
                  }
                }
              }
            }
          }
        }
		
        double sumfetch = 0;
        double sumcos = 0;
        double sumsin =0;
        for(int k = 0; k < 32; k++) {
          double anglesector = pi + k * 2 * pi / 32;
          double cosfetch = cos(anglesector)*minfetch(k)/land1cellszx;
          double sinfetch = sin(anglesector)*minfetch(k)/land1cellszy;
          sumfetch += sqrt((cosfetch*cosfetch) + (sinfetch*sinfetch));
          sumcos += cosfetch;
          sumsin += sinfetch;
          
          //          sumfetch += minfetch(k);
          //          sumcos += cos(anglesector)*minfetch(k);
          //          sumsin += sin(anglesector)*minfetch(k);
        }
        wx32(i,j)=sumfetch;
        //        wx32(i,j+ncolx)=atan2(sumsin,sumcos);
        //        wx32(i,j+ncolx*2)=sqrt(pow((sumsin / sumfetch), 2.0) + pow((sumcos / sumfetch), 2.0));
      }
    }
    if (verbose == 2) Rcout << j << std::endl;
  }
  
  return wx32;
}

//' Compute wave fetch metrics for coastal cells with power-transformed distances (Variant 1b).
//'
//' This function is similar to coastwx32twomapllranf1 but uses a power transformation (dpower)
//' to modify the random distance calculation. It also includes an alternative mapping strategy for cells
//' outside the primary domain.
//'
//' @param land1 IntegerMatrix for the primary land area.
//' @param land2 IntegerMatrix for the secondary land area used for alternative mapping.
//' @param dx Integer specifying the search distance for near-coast identification.
//' @param dwx Integer controlling the magnitude of the random search distance.
//' @param cl Integer indicating the class value for coastal cells.
//' @param nsearches Integer specifying the number of random search iterations.
//' @param dpower Double exponent used for transforming the random distance.
//' @param land1llx Double for the lower left x-coordinate of land1.
//' @param land1uly Double for the upper left y-coordinate of land1.
//' @param land1cellszx Double indicating the cell size in the x-direction for land1.
//' @param land1cellszy Double indicating the cell size in the y-direction for land1.
//' @param land2llx Double for the lower left x-coordinate of land2.
//' @param land2uly Double for the upper left y-coordinate of land2.
//' @param land2cellszx Double indicating the cell size in the x-direction for land2.
//' @param land2cellszy Double indicating the cell size in the y-direction for land2.
//' @param verbose Integer print messages to track progress. Use 0 for no messages, 1 for progress bar or 2 for text messages.
//' @return A NumericMatrix containing the computed wind fetch metric.
// [[Rcpp::export]]
NumericMatrix coastwx32twomapllranf1b(IntegerMatrix land1, IntegerMatrix land2, int dx, int dwx, int cl, int nsearches, double dpower,
                                    double land1llx, double land1uly, double land1cellszx, double land1cellszy, 
                                    double land2llx, double land2uly, double land2cellszx, double land2cellszy,
                                    int verbose = 1) {
  int ncolx = land1.ncol();
  int nrowx = land1.nrow();
  int lncols =land2.ncol();
  int lnrows =land2.nrow();
  
  RProgress::RProgress pb("Processing [:bar] Done: :percent | ETA: :eta", ncolx);
  if (verbose == 2) Rcout << "ncols" << ncolx << "nrows" << nrowx << std::endl;

  double pi = 3.14159265358979323846;
//  int nsearches = 5000;
//  double vdwx = sqrt(dwx)*1000; // changed from int
  double vdwx = pow(dwx,(1/dpower))*1000; // changed from int  pow(10,vdist)
  
  IntegerMatrix coastcell(nrowx,ncolx);
  coastcell = isitnearcoast(land1,dx);
  
  NumericMatrix wx32(nrowx,ncolx);
  
  if (verbose == 1) pb.tick(0);
  for(int j = 0; j < ncolx; j++) {
    if (verbose == 1) pb.tick();
    for(int i = 0; i < nrowx; i++) {
      if (coastcell(i,j)==cl) {
        NumericVector minfetch(32);
        for(int k = 0; k < 32; k++) {
          minfetch(k)=dwx*land1cellszy;
        }
        // Random square root distance and angle
        for(int count = 1; count < nsearches+1;count+=1) {
          
          double vdist = ((double) rand() / RAND_MAX) * vdwx/1000;  // changed from int
          double ranangle = 2 * pi * ((double) rand() / RAND_MAX);
          // double vx = vdist*vdist*sin(vang);
          // double vy = vdist*vdist*cos(vang);
          // 
          // int iranangle = rand() % 10000;
          // int iranlogdist = rand() % 10000;
          
          //double randist = vdist*vdist;
		      double randist = pow(vdist,dpower);
		            
          //Rcout << ranangle << "  " << ranlogdist << "  " << randist << std::endl;
          
          int jshift = randist*sin(ranangle);
          int ishift = randist*cos(ranangle);
          if(j+jshift > -1 && j+jshift < ncolx && i+ishift> -1 && i+ishift<nrowx) {
            if(land1(i+ishift,j+jshift)==1) {
              double yishift=land1cellszy*ishift;
              double xjshift=land1cellszx*jshift;
              double disttoland =  sqrt((xjshift*xjshift)+(yishift*yishift)); 
              int angletoland = 32 * (pi + atan2(xjshift,yishift))/(2*pi); 
              if (angletoland>31) {angletoland -=32;}
              // Rcout << jshift << "  " << ishift << "  " << yishift << "   "<< xjshift << std::endl;
              //  Rcout << ranlogdist << "  " << ranangle << "  " << angletoland << "   "<< disttoland << std::endl;
              if (minfetch(angletoland)>disttoland ) {
                minfetch(angletoland)=disttoland;
              }
            }
          }
          else {
            // alt-map search code here
            int lj = (land1llx - land2llx + land1cellszx * (j + jshift)) / land2cellszx;
            int li = (land2uly - land1uly + land1cellszy * (i + ishift)) / land2cellszy;
            if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
              if(land2(li,lj)==1) {
                double yishift=land1cellszy*ishift;
                double xjshift=land1cellszx*jshift;
                double disttoland = sqrt((xjshift*xjshift)+(yishift*yishift)); 
                int angletoland = 32 * (pi + atan2(xjshift,yishift))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << jshift << "  " << ishift << "  " << yishift << "   "<< xjshift << std::endl;
                // Rcout << ranlogdist << "  " << ranangle << "  " << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
          }
        }
        
		        // No Inner search         
		
        double sumfetch = 0;
        double sumcos = 0;
        double sumsin =0;
        for(int k = 0; k < 32; k++) {
          double anglesector = pi + k * 2 * pi / 32;
          double cosfetch = cos(anglesector)*minfetch(k)/land1cellszx;
          double sinfetch = sin(anglesector)*minfetch(k)/land1cellszy;
          sumfetch += sqrt((cosfetch*cosfetch) + (sinfetch*sinfetch));
          sumcos += cosfetch;
          sumsin += sinfetch;
          
          //          sumfetch += minfetch(k);
          //          sumcos += cos(anglesector)*minfetch(k);
          //          sumsin += sin(anglesector)*minfetch(k);
        }
        wx32(i,j)=sumfetch;
        //        wx32(i,j+ncolx)=atan2(sumsin,sumcos);
        //        wx32(i,j+ncolx*2)=sqrt(pow((sumsin / sumfetch), 2.0) + pow((sumcos / sumfetch), 2.0));
      }
    }
    if (verbose == 2) Rcout << j << std::endl;
  }
  
  return wx32;
}


NumericMatrix coastwx32twomapllranf1a(IntegerMatrix land1, IntegerMatrix land2, int dx, int dwx, int cl, int nsearches, double dpower,
                                    double land1llx, double land1uly, double land1cellszx, double land1cellszy, 
                                    double land2llx, double land2uly, double land2cellszx, double land2cellszy,
                                    int verbose = 1) {
  int ncolx = land1.ncol();
  int nrowx = land1.nrow();
  int lncols =land2.ncol();
  int lnrows =land2.nrow();
  
  if (verbose == 2) Rcout << "ncols" << ncolx << "nrows" << nrowx << std::endl;
  double pi = 3.14159265358979323846;
//  int nsearches = 5000;
//  double vdwx = sqrt(dwx)*1000; // changed from int
  double vdwx = pow(dwx,(1/dpower))*1000; // changed from int  pow(10,vdist)
  
  IntegerMatrix coastcell(nrowx,ncolx);
  coastcell = isitnearcoast(land1,dx);
  
  NumericMatrix wx32(nrowx,ncolx);
  
  for(int j = 0; j < ncolx; j++) {
    for(int i = 0; i < nrowx; i++) {
      if (coastcell(i,j)==cl) {
        NumericVector minfetch(32);
        for(int k = 0; k < 32; k++) {
          minfetch(k)=dwx*land1cellszy;
        }
        // Random square root distance and angle
        for(int count = 1; count < nsearches+1;count+=1) {
          
          double vdist = ((double) rand() / RAND_MAX) * vdwx/1000;  // changed from int
          double ranangle = 2 * pi * ((double) rand() / RAND_MAX);
          // double vx = vdist*vdist*sin(vang);
          // double vy = vdist*vdist*cos(vang);
          // 
          // int iranangle = rand() % 10000;
          // int iranlogdist = rand() % 10000;
          
          //double randist = vdist*vdist;
		  double randist = pow(vdist,dpower);
		            
          //Rcout << ranangle << "  " << ranlogdist << "  " << randist << std::endl;
          
          int jshift = randist*sin(ranangle);
          int ishift = randist*cos(ranangle);
          if(j+jshift > -1 && j+jshift < ncolx && i+ishift> -1 && i+ishift<nrowx) {
            if(land1(i+ishift,j+jshift)==1) {
              double yishift=land1cellszy*ishift;
              double xjshift=land1cellszx*jshift;
              double disttoland =  sqrt((xjshift*xjshift)+(yishift*yishift)); 
              int angletoland = 32 * (pi + atan2(xjshift,yishift))/(2*pi); 
              if (angletoland>31) {angletoland -=32;}
              // Rcout << jshift << "  " << ishift << "  " << yishift << "   "<< xjshift << std::endl;
              //  Rcout << ranlogdist << "  " << ranangle << "  " << angletoland << "   "<< disttoland << std::endl;
              if (minfetch(angletoland)>disttoland ) {
                minfetch(angletoland)=disttoland;
              }
            }
          }
          else {
            // alt-map search code here
            int lj = (land1llx - land2llx + land1cellszx * (j + jshift)) / land2cellszx;
            int li = (land2uly - land1uly + land1cellszy * (i + ishift)) / land2cellszy;
            if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
              if(land2(li,lj)==1) {
                double yishift=land1cellszy*ishift;
                double xjshift=land1cellszx*jshift;
                double disttoland = sqrt((xjshift*xjshift)+(yishift*yishift)); 
                int angletoland = 32 * (pi + atan2(xjshift,yishift))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << jshift << "  " << ishift << "  " << yishift << "   "<< xjshift << std::endl;
                // Rcout << ranlogdist << "  " << ranangle << "  " << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
          }
        }
        
		        // Inner search         
        for(int jshift = -10; jshift < 11; jshift+=1) {
          for(int ishift = -10; ishift < 11; ishift+=1) {
            if(j+jshift > -1 && j+jshift < ncolx && i+ishift> -1 && i+ishift<nrowx) {
              if(land1(i+ishift,j+jshift)==1) {
                double yishiftj=land1cellszy*ishift;
                double xjshiftj=land1cellszx*jshift;
                double disttoland = sqrt((xjshiftj*xjshiftj)+(yishiftj*yishiftj)); 
                int angletoland = 32 * (pi + atan2(xjshiftj,yishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                if (angletoland<0) {angletoland +=32;}
                if (minfetch(angletoland)>disttoland ) {minfetch(angletoland)=disttoland;}
                // for (int ashift = -1; ashift<2; ashift +=1) {
                //   int angleref = angletoland + ashift;
                //   if (angleref>31) {angleref -=32;}
                //   if (angleref<0) {angleref +=32;}
                //   // Rcout << "angleref innerloop" << angleref << "   "<< disttoland << std::endl;
                //   if (minfetch(angleref)>disttoland ) {
                //     minfetch(angleref)=disttoland;
                //   }
                //}
              }
            } 
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellszx * (j + jshift)) / land2cellszx;
              int li = (land2uly - land1uly + land1cellszy * (i + ishift)) / land2cellszy;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double yishiftj=land1cellszy*ishift;
                  double xjshiftj=land1cellszx*jshift;
                  double disttoland = sqrt((xjshiftj*xjshiftj)+(yishiftj*yishiftj)); 
                  int angletoland = 32 * (pi + atan2(xjshiftj,yishiftj))/(2*pi); 
                  if (angletoland>31) {angletoland -=32;}
                  // Rcout << angletoland << "   "<< disttoland << std::endl;
                  if (minfetch(angletoland)>disttoland ) {
                    minfetch(angletoland)=disttoland;
                  }
                }
              }
            }
          }
        }
		
        double sumfetch = 0;
        double sumcos = 0;
        double sumsin =0;
        for(int k = 0; k < 32; k++) {
          double anglesector = pi + k * 2 * pi / 32;
          double cosfetch = cos(anglesector)*minfetch(k)/land1cellszx;
          double sinfetch = sin(anglesector)*minfetch(k)/land1cellszy;
          sumfetch += sqrt((cosfetch*cosfetch) + (sinfetch*sinfetch));
          sumcos += cosfetch;
          sumsin += sinfetch;
          
          //          sumfetch += minfetch(k);
          //          sumcos += cos(anglesector)*minfetch(k);
          //          sumsin += sin(anglesector)*minfetch(k);
        }
        wx32(i,j)=sumfetch;
        //        wx32(i,j+ncolx)=atan2(sumsin,sumcos);
        //        wx32(i,j+ncolx*2)=sqrt(pow((sumsin / sumfetch), 2.0) + pow((sumcos / sumfetch), 2.0));
      }
    }
    if (verbose == 2) Rcout << j << std::endl;
  }
  
  return wx32;
}


//' Compute wind fetch metrics for coastal cells using log-transformed distances (Variant 2).
//'
//' This function calculates a wind fetch metric (wx32) for coastal cells using a log-transformed
//' approach for random distances and includes alternative mapping for cells falling outside the main area.
//'
//' @param land1 IntegerMatrix for the primary land area.
//' @param land2 IntegerMatrix for the secondary land area used for alternative mapping.
//' @param dx Integer specifying the search distance for near-coast identification.
//' @param dwx Integer controlling the magnitude of the random search distance.
//' @param cl Integer indicating the class value for coastal cells.
//' @param nsearches Integer specifying the number of random search iterations.
//' @param land1llx Double for the lower left x-coordinate of land1.
//' @param land1uly Double for the upper left y-coordinate of land1.
//' @param land1cellszx Double indicating the cell size in the x-direction for land1.
//' @param land1cellszy Double indicating the cell size in the y-direction for land1.
//' @param land2llx Double for the lower left x-coordinate of land2.
//' @param land2uly Double for the upper left y-coordinate of land2.
//' @param land2cellszx Double indicating the cell size in the x-direction for land2.
//' @param land2cellszy Double indicating the cell size in the y-direction for land2.
//' @param verbose Integer print messages to track progress. Use 0 for no messages, 1 for progress bar or 2 for text messages.
//' @return A NumericMatrix containing the computed wind fetch metric.
// [[Rcpp::export]]
NumericMatrix coastwx32twomapllranf2(IntegerMatrix land1, IntegerMatrix land2, int dx, int dwx, int cl, int nsearches,
                                     double land1llx, double land1uly, double land1cellszx, double land1cellszy, 
                                     double land2llx, double land2uly, double land2cellszx, double land2cellszy,
                                     int verbose = 1) {
  int ncolx = land1.ncol();
  int nrowx = land1.nrow();
  int lncols =land2.ncol();
  int lnrows =land2.nrow();
  
  RProgress::RProgress pb("Processing [:bar] Done: :percent | ETA: :eta", ncolx);
  if (verbose == 2) Rcout << "ncols" << ncolx << "nrows" << nrowx << std::endl;

  double pi = 3.14159265358979323846;
//  int nsearches = 2000;
  int vdwx = log10(dwx) - 1.0;
  
  IntegerMatrix coastcell(nrowx,ncolx);
  coastcell = isitnearcoast(land1,dx);
  
  NumericMatrix wx32(nrowx,ncolx);
  
  if (verbose == 1) pb.tick(0);
  for(int j = 0; j < ncolx; j++) {
    if (verbose == 1) pb.tick();
    for(int i = 0; i < nrowx; i++) {
      if (coastcell(i,j)==cl) {
        NumericVector minfetch(32);
        for(int k = 0; k < 32; k++) {
          minfetch(k)=dwx*land1cellszy;
        }
        // Random log distance and angle
        for(int count = 1; count < nsearches+1;count+=1) {
          
          double vdist = 1 + vdwx * ((double) rand() / RAND_MAX) ;
          double ranangle = 2 * pi * ((double) rand() / RAND_MAX);
          double randist = pow(10,vdist);
          
          //Rcout << ranangle << "  " << ranlogdist << "  " << randist << std::endl;
          
          int jshift = randist*sin(ranangle);
          int ishift = randist*cos(ranangle);
          if(j+jshift > -1 && j+jshift < ncolx && i+ishift> -1 && i+ishift<nrowx) {
            if(land1(i+ishift,j+jshift)==1) {
              double yishift=land1cellszy*ishift;
              double xjshift=land1cellszx*jshift;
              double disttoland =  sqrt((xjshift*xjshift)+(yishift*yishift)); 
              int angletoland = 32 * (pi + atan2(xjshift,yishift))/(2*pi); 
              if (angletoland>31) {angletoland -=32;}
              // Rcout << jshift << "  " << ishift << "  " << yishift << "   "<< xjshift << std::endl;
              //  Rcout << ranlogdist << "  " << ranangle << "  " << angletoland << "   "<< disttoland << std::endl;
              if (minfetch(angletoland)>disttoland ) {
                minfetch(angletoland)=disttoland;
              }
            }
          }
          else {
            // alt-map search code here
            int lj = (land1llx - land2llx + land1cellszx * (j + jshift)) / land2cellszx;
            int li = (land2uly - land1uly + land1cellszy * (i + ishift)) / land2cellszy;
            if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
              if(land2(li,lj)==1) {
                double yishift=land1cellszy*ishift;
                double xjshift=land1cellszx*jshift;
                double disttoland = sqrt((xjshift*xjshift)+(yishift*yishift)); 
                int angletoland = 32 * (pi + atan2(xjshift,yishift))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << jshift << "  " << ishift << "  " << yishift << "   "<< xjshift << std::endl;
                // Rcout << ranlogdist << "  " << ranangle << "  " << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
          }
        }
        
        // Inner search         
        for(int jshift = -10; jshift < 11; jshift+=1) {
          for(int ishift = -10; ishift < 11; ishift+=1) {
            if(j+jshift > -1 && j+jshift < ncolx && i+ishift> -1 && i+ishift<nrowx) {
              if(land1(i+ishift,j+jshift)==1) {
                double yishiftj=land1cellszy*ishift;
                double xjshiftj=land1cellszx*jshift;
                double disttoland = sqrt((xjshiftj*xjshiftj)+(yishiftj*yishiftj)); 
                int angletoland = 32 * (pi + atan2(xjshiftj,yishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                if (angletoland<0) {angletoland +=32;}
                if (minfetch(angletoland)>disttoland ) {minfetch(angletoland)=disttoland;}
                // for (int ashift = -1; ashift<2; ashift +=1) {
                //   int angleref = angletoland + ashift;
                //   if (angleref>31) {angleref -=32;}
                //   if (angleref<0) {angleref +=32;}
                //   // Rcout << "angleref innerloop" << angleref << "   "<< disttoland << std::endl;
                //   if (minfetch(angleref)>disttoland ) {
                //     minfetch(angleref)=disttoland;
                //   }
                //}
              }
            } 
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellszx * (j + jshift)) / land2cellszx;
              int li = (land2uly - land1uly + land1cellszy * (i + ishift)) / land2cellszy;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double yishiftj=land1cellszy*ishift;
                  double xjshiftj=land1cellszx*jshift;
                  double disttoland = sqrt((xjshiftj*xjshiftj)+(yishiftj*yishiftj)); 
                  int angletoland = 32 * (pi + atan2(xjshiftj,yishiftj))/(2*pi); 
                  if (angletoland>31) {angletoland -=32;}
                  // Rcout << angletoland << "   "<< disttoland << std::endl;
                  if (minfetch(angletoland)>disttoland ) {
                    minfetch(angletoland)=disttoland;
                  }
                }
              }
            }
          }
        }       
        
        
        double sumfetch = 0;
        double sumcos = 0;
        double sumsin =0;
        for(int k = 0; k < 32; k++) {
          double anglesector = pi + k * 2 * pi / 32;
          double cosfetch = cos(anglesector)*minfetch(k)/land1cellszx;
          double sinfetch = sin(anglesector)*minfetch(k)/land1cellszy;
          sumfetch += sqrt((cosfetch*cosfetch) + (sinfetch*sinfetch));
          sumcos += cosfetch;
          sumsin += sinfetch;
          
          //          sumfetch += minfetch(k);
          //          sumcos += cos(anglesector)*minfetch(k);
          //          sumsin += sin(anglesector)*minfetch(k);
        }
        wx32(i,j)=sumfetch;
        //        wx32(i,j+ncolx)=atan2(sumsin,sumcos);
        //        wx32(i,j+ncolx*2)=sqrt(pow((sumsin / sumfetch), 2.0) + pow((sumcos / sumfetch), 2.0));
      }
    }
    if (verbose == 2) Rcout << j << std::endl;
  }
  
  return wx32;
}

