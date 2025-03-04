/* Wave fetch models in C++ for R */

#include <Rcpp.h>
#include <cmath>
#include <RProgress.h>
// [[Rcpp::depends("progress")]]
using namespace Rcpp;



// Compute the mean of a numeric vector.
//
// @param x NumericVector of values.
// @return The arithmetic mean of x.
// [[Rcpp::export]]
double meanC(NumericVector x) {
  int n = x.size();
  double total = 0;
  
  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total / n;
}



// Identify coastal cells in a land matrix.
//
// This function examines each cell and its 3×3 neighborhood in the provided land matrix.
// A cell is considered coastal if it is land (value 1) but not completely surrounded by land.
//
// @param land IntegerMatrix representing land (1) and water (0).
// @return An IntegerMatrix with coastal cells marked as 1.
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



// Identify near-coast cells based on a land matrix and a search distance.
//
// For cells identified as coastal, this function marks surrounding water cells (value 0) 
// within a distance dx as near-coast (value 2).
//
// @param land IntegerMatrix representing land (1) and water (0).
// @param dx Integer specifying the search distance (in cells).
// @return An IntegerMatrix with coastal cells marked as 1 and near-coast cells marked as 2.
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


//' Compute the wave fetch metric using log-transformed distances.
//'
//' This function calculates a wave fetch metric (wx32) for coastal cells,
//' and uses a log-transformed approach for the random distance search.
//' An alternative mapping search is performed when the primary search is out-of-bound.
//'
//' The original name of this function was `coastwx32twomaplatlonf`
//'
//' @param land1 IntegerMatrix representing the primary land area.
//' @param land2 IntegerMatrix representing the secondary land area used for alternative mapping.
//' @param dx Integer specifying the search distance for near-coast detection.
//' @param dwx Integer controlling the magnitude of the random search distance.
//' @param cl Integer indicating the coastal cell class to process.
//' @param nsearches Integer specifying the number of iterations for the random search.
//' @param land1llx Double: lower left x-coordinate for land1.
//' @param land1uly Double: upper left y-coordinate for land1.
//' @param land1cellszx Double: cell size in the x-direction for land1.
//' @param land1cellszy Double: cell size in the y-direction for land1.
//' @param land2llx Double: lower left x-coordinate for land2.
//' @param land2uly Double: upper left y-coordinate for land2.
//' @param land2cellszx Double: cell size in the x-direction for land2.
//' @param land2cellszy Double: cell size in the y-direction for land2.
//' @param verbose Integer print messages to track progress. Use 0 for no messages, 1 for progress bar or 2 for text messages.
//' @return A NumericMatrix with the computed wave fetch metric based on the log-transformed search approach.
//' @export
// [[Rcpp::export]]
NumericMatrix coastal_wave_fetch(IntegerMatrix land1, IntegerMatrix land2, int dx, int dwx, int cl, int jit,
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
  
  if (verbose >= 1) Rcout << "Identifying near coast cells\n\n" << std::endl;
  IntegerMatrix coastcell(nrowx,ncolx);
  coastcell = isitnearcoast(land1,dx);
  
  NumericMatrix wx32(nrowx,ncolx);
  
  if (verbose >= 1) Rcout << "Calculating metric...\n\n" << std::endl;
  if (verbose == 1) pb.tick(0);
  for(int j = 0; j < ncolx; j++) {
    if (verbose == 1) pb.tick();
    for(int i = 0; i < nrowx; i++) {
      if (coastcell(i,j)==cl) {
        NumericVector minfetch(32);
        for(int k = 0; k < 32; k++) {
          minfetch(k)=dwx*land1cellszy;
        }
        // Outer search        
        for(int jshift = -dwx; jshift < dwx+1; jshift+=100) {
          for(int ishift = -dwx; ishift < dwx+1; ishift+=100) {
            int jshiftj = jshift + jit*(-50 + rand() % 100);
            int ishiftj = ishift + jit*(-50 + rand() % 100);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land1(i+ishiftj,j+jshiftj)==1) {
                double yishiftj=land1cellszy*ishiftj;
                double xjshiftj=land1cellszx*jshiftj;
                double disttoland = sqrt((xjshiftj*xjshiftj)+(yishiftj*yishiftj)); 
                int angletoland = 32 * (pi + atan2(xjshiftj,yishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellszx * (j + jshiftj)) / land2cellszx;
              int li = (land2uly - land1uly + land1cellszy * (i + ishiftj)) / land2cellszy;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double yishiftj=land1cellszy*ishiftj;
                  double xjshiftj=land1cellszx*jshiftj;
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
        // Middle search
        for(int jshift = -0.1*dwx; jshift < 0.1*dwx+1; jshift+=10) {
          for(int ishift = -0.1*dwx; ishift < 0.1*dwx+1; ishift+=10) {
            int jshiftj = jshift + jit*(-5 + rand() % 10);
            int ishiftj = ishift + jit*(-5 + rand() % 10);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land1(i+ishiftj,j+jshiftj)==1) {
                double yishiftj=land1cellszy*ishiftj;
                double xjshiftj=land1cellszx*jshiftj;
                double disttoland = sqrt((xjshiftj*xjshiftj)+(yishiftj*yishiftj)); 
                int angletoland = 32 * (pi + atan2(xjshiftj,yishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellszx * (j + jshiftj)) / land2cellszx;
              int li = (land2uly - land1uly + land1cellszy * (i + ishiftj)) / land2cellszy;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double yishiftj=land1cellszy*ishiftj;
                  double xjshiftj=land1cellszx*jshiftj;
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
                // Rcout << angletoland << "   "<< disttoland << std::endl;
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
          
        }
        wx32(i,j)=sumfetch;
        //        wx32(i,j+ncolx)=atan2(sumsin,sumcos);
        //        wx32(i,j+ncolx*2)=sqrt(pow((sumsin / sumfetch), 2.0) + pow((sumcos / sumfetch), 2.0));
      }
    }
    if (verbose == 2) Rcout << j << std::endl;
  }
  if (verbose >= 1) Rcout << "Done!" << std::endl;
  return wx32;
}



//' Compute the wave orientation metric using log-transformed distances.
//'
//' This function calculates a wave orientation metric (wx32) for coastal cells,
//' and uses a log-transformed approach for the random distance search.
//' An alternative mapping search is performed when the primary search is out-of-bound.
//'
//' The original name of this function was `coastwx32twomaplatlono`
//'
//' @param land1 IntegerMatrix representing the primary land area.
//' @param land2 IntegerMatrix representing the secondary land area used for alternative mapping.
//' @param dx Integer specifying the search distance for near-coast detection.
//' @param dwx Integer controlling the magnitude of the random search distance.
//' @param cl Integer indicating the coastal cell class to process.
//' @param nsearches Integer specifying the number of iterations for the random search.
//' @param land1llx Double: lower left x-coordinate for land1.
//' @param land1uly Double: upper left y-coordinate for land1.
//' @param land1cellszx Double: cell size in the x-direction for land1.
//' @param land1cellszy Double: cell size in the y-direction for land1.
//' @param land2llx Double: lower left x-coordinate for land2.
//' @param land2uly Double: upper left y-coordinate for land2.
//' @param land2cellszx Double: cell size in the x-direction for land2.
//' @param land2cellszy Double: cell size in the y-direction for land2.
//' @param verbose Integer print messages to track progress. Use 0 for no messages, 1 for progress bar or 2 for text messages.
//' @return A NumericMatrix with the computed wave fetch metric based on the log-transformed search approach.
//' @export
// [[Rcpp::export]]
NumericMatrix coastal_wave_orientation(IntegerMatrix land1, IntegerMatrix land2, int dx, int dwx, int cl, int jit,
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
  
  if (verbose >= 1) Rcout << "Identifying near coast cells\n\n" << std::endl;
  IntegerMatrix coastcell(nrowx,ncolx);
  coastcell = isitnearcoast(land1,dx);
  
  NumericMatrix wx32(nrowx,ncolx);
  
  if (verbose >= 1) Rcout << "Calculating metric...\n\n" << std::endl;
  if (verbose == 1) pb.tick(0);
  for(int j = 0; j < ncolx; j++) {
    if (verbose == 1) pb.tick();
    for(int i = 0; i < nrowx; i++) {
      if (coastcell(i,j)==cl) {
        NumericVector minfetch(32);
        for(int k = 0; k < 32; k++) {
          minfetch(k)=dwx*land1cellszy;
        }
        // Outer search        
        for(int jshift = -dwx; jshift < dwx+1; jshift+=100) {
          for(int ishift = -dwx; ishift < dwx+1; ishift+=100) {
            int jshiftj = jshift + jit*(-50 + rand() % 100);
            int ishiftj = ishift + jit*(-50 + rand() % 100);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land1(i+ishiftj,j+jshiftj)==1) {
                double yishiftj=land1cellszy*ishiftj;
                double xjshiftj=land1cellszx*jshiftj;
                double disttoland = sqrt((xjshiftj*xjshiftj)+(yishiftj*yishiftj)); 
                int angletoland = 32 * (pi + atan2(xjshiftj,yishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellszx * (j + jshiftj)) / land2cellszx;
              int li = (land2uly - land1uly + land1cellszy * (i + ishiftj)) / land2cellszy;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double yishiftj=land1cellszy*ishiftj;
                  double xjshiftj=land1cellszx*jshiftj;
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
        // Middle search
        for(int jshift = -0.1*dwx; jshift < 0.1*dwx+1; jshift+=10) {
          for(int ishift = -0.1*dwx; ishift < 0.1*dwx+1; ishift+=10) {
            int jshiftj = jshift + jit*(-5 + rand() % 10);
            int ishiftj = ishift + jit*(-5 + rand() % 10);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land1(i+ishiftj,j+jshiftj)==1) {
                double yishiftj=land1cellszy*ishiftj;
                double xjshiftj=land1cellszx*jshiftj;
                double disttoland = sqrt((xjshiftj*xjshiftj)+(yishiftj*yishiftj)); 
                int angletoland = 32 * (pi + atan2(xjshiftj,yishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellszx * (j + jshiftj)) / land2cellszx;
              int li = (land2uly - land1uly + land1cellszy * (i + ishiftj)) / land2cellszy;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double yishiftj=land1cellszy*ishiftj;
                  double xjshiftj=land1cellszx*jshiftj;
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
                // Rcout << angletoland << "   "<< disttoland << std::endl;
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
          
        }
        //      wx32(i,j)=sumfetch;
        wx32(i,j)=atan2(sumsin,sumcos);
        //        wx32(i,j+ncolx*2)=sqrt(pow((sumsin / sumfetch), 2.0) + pow((sumcos / sumfetch), 2.0));
      }
    }
    if (verbose == 2) Rcout << j << std::endl;
  }
  if (verbose >= 1) Rcout << "Done!" << std::endl;
  return wx32;
}


//' Compute the wave directionality metric using log-transformed distances.
//'
//' This function calculates a wave directionality metric (wx32) for coastal cells,
//' and uses a log-transformed approach for the random distance search.
//' An alternative mapping search is performed when the primary search is out-of-bound.
//'
//' The original name of this function was `coastwx32twomaplatlonv`
//'
//' @param land1 IntegerMatrix representing the primary land area.
//' @param land2 IntegerMatrix representing the secondary land area used for alternative mapping.
//' @param dx Integer specifying the search distance for near-coast detection.
//' @param dwx Integer controlling the magnitude of the random search distance.
//' @param cl Integer indicating the coastal cell class to process.
//' @param nsearches Integer specifying the number of iterations for the random search.
//' @param land1llx Double: lower left x-coordinate for land1.
//' @param land1uly Double: upper left y-coordinate for land1.
//' @param land1cellszx Double: cell size in the x-direction for land1.
//' @param land1cellszy Double: cell size in the y-direction for land1.
//' @param land2llx Double: lower left x-coordinate for land2.
//' @param land2uly Double: upper left y-coordinate for land2.
//' @param land2cellszx Double: cell size in the x-direction for land2.
//' @param land2cellszy Double: cell size in the y-direction for land2.
//' @param verbose Integer print messages to track progress. Use 0 for no messages, 1 for progress bar or 2 for text messages.
//' @return A NumericMatrix with the computed wave fetch metric based on the log-transformed search approach.
//' @export
// [[Rcpp::export]]
NumericMatrix coastal_wave_direction(IntegerMatrix land1, IntegerMatrix land2, int dx, int dwx, int cl, int jit,
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
  
  if (verbose >= 1) Rcout << "Identifying near coast cells\n\n" << std::endl;
  IntegerMatrix coastcell(nrowx,ncolx);
  coastcell = isitnearcoast(land1,dx);
  
  NumericMatrix wx32(nrowx,ncolx);
  
  if (verbose >= 1) Rcout << "Calculating metric...\n\n" << std::endl;
  if (verbose == 1) pb.tick(0);
  for(int j = 0; j < ncolx; j++) {
    if (verbose == 1) pb.tick();
    for(int i = 0; i < nrowx; i++) {
      if (coastcell(i,j)==cl) {
        NumericVector minfetch(32);
        for(int k = 0; k < 32; k++) {
          minfetch(k)=dwx*land1cellszy;
        }
        // Outer search        
        for(int jshift = -dwx; jshift < dwx+1; jshift+=100) {
          for(int ishift = -dwx; ishift < dwx+1; ishift+=100) {
            int jshiftj = jshift + jit*(-50 + rand() % 100);
            int ishiftj = ishift + jit*(-50 + rand() % 100);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land1(i+ishiftj,j+jshiftj)==1) {
                double yishiftj=land1cellszy*ishiftj;
                double xjshiftj=land1cellszx*jshiftj;
                double disttoland = sqrt((xjshiftj*xjshiftj)+(yishiftj*yishiftj)); 
                int angletoland = 32 * (pi + atan2(xjshiftj,yishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellszx * (j + jshiftj)) / land2cellszx;
              int li = (land2uly - land1uly + land1cellszy * (i + ishiftj)) / land2cellszy;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double yishiftj=land1cellszy*ishiftj;
                  double xjshiftj=land1cellszx*jshiftj;
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
        // Middle search
        for(int jshift = -0.1*dwx; jshift < 0.1*dwx+1; jshift+=10) {
          for(int ishift = -0.1*dwx; ishift < 0.1*dwx+1; ishift+=10) {
            int jshiftj = jshift + jit*(-5 + rand() % 10);
            int ishiftj = ishift + jit*(-5 + rand() % 10);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land1(i+ishiftj,j+jshiftj)==1) {
                double yishiftj=land1cellszy*ishiftj;
                double xjshiftj=land1cellszx*jshiftj;
                double disttoland = sqrt((xjshiftj*xjshiftj)+(yishiftj*yishiftj)); 
                int angletoland = 32 * (pi + atan2(xjshiftj,yishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellszx * (j + jshiftj)) / land2cellszx;
              int li = (land2uly - land1uly + land1cellszy * (i + ishiftj)) / land2cellszy;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double yishiftj=land1cellszy*ishiftj;
                  double xjshiftj=land1cellszx*jshiftj;
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
                // Rcout << angletoland << "   "<< disttoland << std::endl;
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
        }
        //      wx32(i,j)=sumfetch;
        //wx32(i,j)=atan2(sumsin,sumcos);
        wx32(i,j)=sqrt(pow((sumsin / sumfetch), 2.0) + pow((sumcos / sumfetch), 2.0));
      }
    }
    if(verbose == 2) Rcout << j << std::endl;
  }
  if (verbose >= 1) Rcout << "Done!" << std::endl;
  return wx32;
}








/* OLD CODE */
NumericMatrix coastwx32(IntegerMatrix land, int dx, int dwx, int cl, int jit) {
  int ncolx = land.ncol();
  int nrowx = land.nrow();
  Rcout << "ncols" << ncolx << "nrows" << nrowx << std::endl;
  double pi = 3.14159265358979323846;
  
  IntegerMatrix coastcell(nrowx,ncolx);
  coastcell = isitnearcoast(land,dx);
  
  NumericMatrix wx32(nrowx,ncolx);
  
  for(int j = 0; j < ncolx; j++) {
    for(int i = 0; i < nrowx; i++) {
      if (coastcell(i,j)==cl) {
        NumericVector minfetch(32);
        for(int k = 0; k < 32; k++) {
          minfetch(k)=dwx;
        }
        // Outer search        
        for(int jshift = -dwx; jshift < dwx+1; jshift+=100) {
          for(int ishift = -dwx; ishift < dwx+1; ishift+=100) {
            int jshiftj = jshift + jit*(-50 + rand() % 100);
            int ishiftj = ishift + jit*(-50 + rand() % 100);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land(i+ishiftj,j+jshiftj)==1) {
                double disttoland = sqrt((jshiftj*jshiftj)+(ishiftj*ishiftj)); 
                int angletoland = 32 * (pi + atan2(jshiftj,ishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
          }
        }
        // Middle search
        for(int jshift = -0.1*dwx; jshift < 0.1*dwx+1; jshift+=10) {
          for(int ishift = -0.1*dwx; ishift < 0.1*dwx+1; ishift+=10) {
            int jshiftj = jshift + jit*(-5 + rand() % 10);
            int ishiftj = ishift + jit*(-5 + rand() % 10);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land(i+ishiftj,j+jshiftj)==1) {
                double disttoland = sqrt((jshiftj*jshiftj)+(ishiftj*ishiftj)); 
                int angletoland = 32 * (pi + atan2(jshiftj,ishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
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
              if(land(i+ishift,j+jshift)==1) {
                double disttoland = sqrt((jshift*jshift)+(ishift*ishift)); 
                int angletoland = 32 * (pi + atan2(jshift,ishift))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
          }
        }       
        double sumfetch = 0;
        for(int k = 0; k < 32; k++) {
          sumfetch += minfetch(k);
        }
        wx32(i,j)=sumfetch;
      }
    }
    Rcout << j << std::endl;
  }
  return wx32;
}

NumericMatrix coastwx32mod(IntegerMatrix land, int dx, int dwx, int outstep, int midstep, int cl, int jit) {
  int ncolx = land.ncol();
  int nrowx = land.nrow();
  Rcout << "ncols" << ncolx << "nrows" << nrowx << std::endl;
  double pi = 3.14159265358979323846;
  
  IntegerMatrix coastcell(nrowx,ncolx);
  coastcell = isitnearcoast(land,dx);
  
  NumericMatrix wx32(nrowx,ncolx);
  
  for(int j = 0; j < ncolx; j++) {
    for(int i = 0; i < nrowx; i++) {
      if (coastcell(i,j)==cl) {
        NumericVector minfetch(32);
        for(int k = 0; k < 32; k++) {
          minfetch(k)=dwx;
        }
        // Outer search        
        for(int jshift = -dwx; jshift < dwx+1; jshift+=outstep) {
          for(int ishift = -dwx; ishift < dwx+1; ishift+=outstep) {
            int jshiftj = jshift + jit*(-outstep*0.5 + rand() % outstep);
            int ishiftj = ishift + jit*(-outstep*0.5 + rand() % outstep);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land(i+ishiftj,j+jshiftj)==1) {
                double disttoland = sqrt((jshiftj*jshiftj)+(ishiftj*ishiftj)); 
                int angletoland = 32 * (pi + atan2(jshiftj,ishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
          }
        }
        // Middle search
        for(int jshift = -0.1*dwx; jshift < 0.1*dwx+1; jshift+=midstep) {
          for(int ishift = -0.1*dwx; ishift < 0.1*dwx+1; ishift+=midstep) {
            int jshiftj = jshift + jit*(-0.5*midstep + rand() % midstep);
            int ishiftj = ishift + jit*(-0.5*midstep + rand() % midstep);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land(i+ishiftj,j+jshiftj)==1) {
                double disttoland = sqrt((jshiftj*jshiftj)+(ishiftj*ishiftj)); 
                int angletoland = 32 * (pi + atan2(jshiftj,ishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
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
              if(land(i+ishift,j+jshift)==1) {
                double disttoland = sqrt((jshift*jshift)+(ishift*ishift)); 
                int angletoland = 32 * (pi + atan2(jshift,ishift))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
          }
        }       
        double sumfetch = 0;
        for(int k = 0; k < 32; k++) {
          sumfetch += minfetch(k);
        }
        wx32(i,j)=sumfetch;
      }
    }
    Rcout << j << std::endl;
  }
  return wx32;
}


NumericMatrix coastwx32twomap(IntegerMatrix land1, IntegerMatrix land2, int dx, int dwx, int cl, int jit,
                              double land1llx, double land1uly, double land1cellsz, 
                              double land2llx, double land2uly, double land2cellsz) {
  int ncolx = land1.ncol();
  int nrowx = land1.nrow();
  int lncols =land2.ncol();
  int lnrows =land2.nrow();
  
  Rcout << "ncols" << ncolx << "nrows" << nrowx << std::endl;
  double pi = 3.14159265358979323846;
  
  IntegerMatrix coastcell(nrowx,ncolx);
  coastcell = isitnearcoast(land1,dx);
  
  NumericMatrix wx32(nrowx,ncolx);
  
  for(int j = 0; j < ncolx; j++) {
    for(int i = 0; i < nrowx; i++) {
      if (coastcell(i,j)==cl) {
        NumericVector minfetch(32);
        for(int k = 0; k < 32; k++) {
          minfetch(k)=dwx;
        }
        // Outer search        
        for(int jshift = -dwx; jshift < dwx+1; jshift+=100) {
          for(int ishift = -dwx; ishift < dwx+1; ishift+=100) {
            int jshiftj = jshift + jit*(-50 + rand() % 100);
            int ishiftj = ishift + jit*(-50 + rand() % 100);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land1(i+ishiftj,j+jshiftj)==1) {
                double disttoland = sqrt((jshiftj*jshiftj)+(ishiftj*ishiftj)); 
                int angletoland = 32 * (pi + atan2(jshiftj,ishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellsz * (j + jshiftj)) / land2cellsz;
              int li = (land2uly - land1uly + land1cellsz * (i + ishiftj)) / land2cellsz;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double disttoland = sqrt((jshiftj*jshiftj)+(ishiftj*ishiftj)); 
                  int angletoland = 32 * (pi + atan2(jshiftj,ishiftj))/(2*pi); 
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
        // Middle search
        for(int jshift = -0.1*dwx; jshift < 0.1*dwx+1; jshift+=10) {
          for(int ishift = -0.1*dwx; ishift < 0.1*dwx+1; ishift+=10) {
            int jshiftj = jshift + jit*(-5 + rand() % 10);
            int ishiftj = ishift + jit*(-5 + rand() % 10);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land1(i+ishiftj,j+jshiftj)==1) {
                double disttoland = sqrt((jshiftj*jshiftj)+(ishiftj*ishiftj)); 
                int angletoland = 32 * (pi + atan2(jshiftj,ishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellsz * (j + jshiftj)) / land2cellsz;
              int li = (land2uly - land1uly + land1cellsz * (i + ishiftj)) / land2cellsz;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double disttoland = sqrt((jshiftj*jshiftj)+(ishiftj*ishiftj)); 
                  int angletoland = 32 * (pi + atan2(jshiftj,ishiftj))/(2*pi); 
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
        
        // Inner search         
        for(int jshift = -10; jshift < 11; jshift+=1) {
          for(int ishift = -10; ishift < 11; ishift+=1) {
            if(j+jshift > -1 && j+jshift < ncolx && i+ishift> -1 && i+ishift<nrowx) {
              if(land1(i+ishift,j+jshift)==1) {
                double disttoland = sqrt((jshift*jshift)+(ishift*ishift)); 
                int angletoland = 32 * (pi + atan2(jshift,ishift))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            } 
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellsz * (j + jshift)) / land2cellsz;
              int li = (land2uly - land1uly + land1cellsz * (i + ishift)) / land2cellsz;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double disttoland = sqrt((jshift*jshift)+(ishift*ishift)); 
                  int angletoland = 32 * (pi + atan2(jshift,ishift))/(2*pi); 
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
        for(int k = 0; k < 32; k++) {
          sumfetch += minfetch(k);
        }
        wx32(i,j)=sumfetch;
      }
    }
    Rcout << j << std::endl;
  }
  return wx32;
}



NumericMatrix coastwx32twomaplatlon(IntegerMatrix land1, IntegerMatrix land2, int dx, int dwx, int cl, int jit,
                                    double land1llx, double land1uly, double land1cellszx, double land1cellszy, 
                                    double land2llx, double land2uly, double land2cellszx, double land2cellszy) {
  int ncolx = land1.ncol();
  int nrowx = land1.nrow();
  int lncols =land2.ncol();
  int lnrows =land2.nrow();
  
  Rcout << "ncols" << ncolx << "nrows" << nrowx << std::endl;
  double pi = 3.14159265358979323846;
  
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
        // Outer search        
        for(int jshift = -dwx; jshift < dwx+1; jshift+=100) {
          for(int ishift = -dwx; ishift < dwx+1; ishift+=100) {
            int jshiftj = jshift + jit*(-50 + rand() % 100);
            int ishiftj = ishift + jit*(-50 + rand() % 100);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land1(i+ishiftj,j+jshiftj)==1) {
                double yishiftj=land1cellszy*ishiftj;
                double xjshiftj=land1cellszx*jshiftj;
                double disttoland = sqrt((xjshiftj*xjshiftj)+(yishiftj*yishiftj)); 
                int angletoland = 32 * (pi + atan2(xjshiftj,yishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellszx * (j + jshiftj)) / land2cellszx;
              int li = (land2uly - land1uly + land1cellszy * (i + ishiftj)) / land2cellszy;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double yishiftj=land1cellszy*ishiftj;
                  double xjshiftj=land1cellszx*jshiftj;
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
        // Middle search
        for(int jshift = -0.1*dwx; jshift < 0.1*dwx+1; jshift+=10) {
          for(int ishift = -0.1*dwx; ishift < 0.1*dwx+1; ishift+=10) {
            int jshiftj = jshift + jit*(-5 + rand() % 10);
            int ishiftj = ishift + jit*(-5 + rand() % 10);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land1(i+ishiftj,j+jshiftj)==1) {
                double yishiftj=land1cellszy*ishiftj;
                double xjshiftj=land1cellszx*jshiftj;
                double disttoland = sqrt((xjshiftj*xjshiftj)+(yishiftj*yishiftj)); 
                int angletoland = 32 * (pi + atan2(xjshiftj,yishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellszx * (j + jshiftj)) / land2cellszx;
              int li = (land2uly - land1uly + land1cellszy * (i + ishiftj)) / land2cellszy;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double yishiftj=land1cellszy*ishiftj;
                  double xjshiftj=land1cellszx*jshiftj;
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
                // Rcout << angletoland << "   "<< disttoland << std::endl;
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
        for(int k = 0; k < 32; k++) {
          sumfetch += minfetch(k);
        }
        wx32(i,j)=sumfetch;
      }
    }
    Rcout << j << std::endl;
  }
  return wx32;
}


NumericMatrix coastwx32twomap3(IntegerMatrix land1, IntegerMatrix land2, int dx, int dwx, int cl, int jit,
                               double land1llx, double land1uly, double land1cellsz, 
                               double land2llx, double land2uly, double land2cellsz) {
  int ncolx = land1.ncol();
  int nrowx = land1.nrow();
  int lncols =land2.ncol();
  int lnrows =land2.nrow();
  
  Rcout << "ncols" << ncolx << "nrows" << nrowx << std::endl;
  double pi = 3.14159265358979323846;
  
  IntegerMatrix coastcell(nrowx,ncolx);
  coastcell = isitnearcoast(land1,dx);
  
  NumericMatrix wx32(nrowx,ncolx*3);
  
  for(int j = 0; j < ncolx; j++) {
    for(int i = 0; i < nrowx; i++) {
      if (coastcell(i,j)==cl) {
        NumericVector minfetch(32);
        for(int k = 0; k < 32; k++) {
          minfetch(k)=dwx;
        }
        // Outer search        
        for(int jshift = -dwx; jshift < dwx+1; jshift+=100) {
          for(int ishift = -dwx; ishift < dwx+1; ishift+=100) {
            int jshiftj = jshift + jit*(-50 + rand() % 100);
            int ishiftj = ishift + jit*(-50 + rand() % 100);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land1(i+ishiftj,j+jshiftj)==1) {
                double disttoland = sqrt((jshiftj*jshiftj)+(ishiftj*ishiftj)); 
                int angletoland = 32 * (pi + atan2(jshiftj,ishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellsz * (j + jshiftj)) / land2cellsz;
              int li = (land2uly - land1uly + land1cellsz * (i + ishiftj)) / land2cellsz;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double disttoland = sqrt((jshiftj*jshiftj)+(ishiftj*ishiftj)); 
                  int angletoland = 32 * (pi + atan2(jshiftj,ishiftj))/(2*pi); 
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
        // Middle search
        for(int jshift = -0.1*dwx; jshift < 0.1*dwx+1; jshift+=6) {
          for(int ishift = -0.1*dwx; ishift < 0.1*dwx+1; ishift+=6) {
            int jshiftj = jshift + jit*(-3 + rand() % 6);
            int ishiftj = ishift + jit*(-3 + rand() % 6);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land1(i+ishiftj,j+jshiftj)==1) {
                double disttoland = sqrt((jshiftj*jshiftj)+(ishiftj*ishiftj)); 
                int angletoland = 32 * (pi + atan2(jshiftj,ishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellsz * (j + jshiftj)) / land2cellsz;
              int li = (land2uly - land1uly + land1cellsz * (i + ishiftj)) / land2cellsz;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double disttoland = sqrt((jshiftj*jshiftj)+(ishiftj*ishiftj)); 
                  int angletoland = 32 * (pi + atan2(jshiftj,ishiftj))/(2*pi); 
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
        
        // Inner search         
        for(int jshift = -10; jshift < 11; jshift+=1) {
          for(int ishift = -10; ishift < 11; ishift+=1) {
            if(j+jshift > -1 && j+jshift < ncolx && i+ishift> -1 && i+ishift<nrowx) {
              if(land1(i+ishift,j+jshift)==1) {
                double disttoland = sqrt((jshift*jshift)+(ishift*ishift)); 
                int angletoland = 32 * (pi + atan2(jshift,ishift))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            } 
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellsz * (j + jshift)) / land2cellsz;
              int li = (land2uly - land1uly + land1cellsz * (i + ishift)) / land2cellsz;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double disttoland = sqrt((jshift*jshift)+(ishift*ishift)); 
                  int angletoland = 32 * (pi + atan2(jshift,ishift))/(2*pi); 
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
          sumfetch += minfetch(k);
          sumcos += cos(anglesector)*minfetch(k);
          sumsin += sin(anglesector)*minfetch(k);
        }
        wx32(i,j)=sumfetch;
        wx32(i,j+ncolx)=atan2(sumsin,sumcos);
        wx32(i,j+ncolx*2)=sqrt(pow((sumsin / sumfetch), 2.0) + pow((sumcos / sumfetch), 2.0));
      }
    }
    Rcout << j << std::endl;
  }
  return wx32;
}


NumericMatrix coastwx32twomap3q(IntegerMatrix land1, IntegerMatrix land2, int dx, int dwx, int cl, int jit,
                                double land1llx, double land1uly, double land1cellsz, 
                                double land2llx, double land2uly, double land2cellsz) {
  int ncolx = land1.ncol();
  int nrowx = land1.nrow();
  int lncols =land2.ncol();
  int lnrows =land2.nrow();
  
  Rcout << "ncols" << ncolx << "nrows" << nrowx << std::endl;
  double pi = 3.14159265358979323846;
  
  IntegerMatrix coastcell(nrowx,ncolx);
  coastcell = isitnearcoast(land1,dx);
  
  NumericMatrix wx32(nrowx,ncolx*3);
  int outstep = ncolx / 20;
  int midstep = ncolx / 200;
  
  for(int j = 0; j < ncolx; j++) {
    for(int i = 0; i < nrowx; i++) {
      if (coastcell(i,j)==cl) {
        NumericVector minfetch(32);
        for(int k = 0; k < 32; k++) {
          minfetch(k)=dwx;
        }
        // Outer search        
        for(int jshift = -dwx; jshift < dwx+1; jshift+=outstep) {
          for(int ishift = -dwx; ishift < dwx+1; ishift+=outstep) {
            int jshiftj = jshift + jit*(-0.5*outstep + rand() % outstep);
            int ishiftj = ishift + jit*(-0.5*outstep + rand() % outstep);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land1(i+ishiftj,j+jshiftj)==1) {
                double disttoland = sqrt((jshiftj*jshiftj)+(ishiftj*ishiftj)); 
                int angletoland = 32 * (pi + atan2(jshiftj,ishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellsz * (j + jshiftj)) / land2cellsz;
              int li = (land2uly - land1uly + land1cellsz * (i + ishiftj)) / land2cellsz;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double disttoland = sqrt((jshiftj*jshiftj)+(ishiftj*ishiftj)); 
                  int angletoland = 32 * (pi + atan2(jshiftj,ishiftj))/(2*pi); 
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
        // Middle search
        for(int jshift = -0.1*dwx; jshift < 0.1*dwx+1; jshift+=midstep) {
          for(int ishift = -0.1*dwx; ishift < 0.1*dwx+1; ishift+=midstep) {
            // for(int jshift = -outstep; jshift < outstep+1; jshift+=midstep) {
            //   for(int ishift = -outstep; ishift < outstep+1; ishift+=midstep) {
            
            int jshiftj = jshift + jit*(-0.5*midstep + rand() % midstep);
            int ishiftj = ishift + jit*(-0.5*midstep + rand() % midstep);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land1(i+ishiftj,j+jshiftj)==1) {
                double disttoland = sqrt((jshiftj*jshiftj)+(ishiftj*ishiftj)); 
                int angletoland = 32 * (pi + atan2(jshiftj,ishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellsz * (j + jshiftj)) / land2cellsz;
              int li = (land2uly - land1uly + land1cellsz * (i + ishiftj)) / land2cellsz;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double disttoland = sqrt((jshiftj*jshiftj)+(ishiftj*ishiftj)); 
                  int angletoland = 32 * (pi + atan2(jshiftj,ishiftj))/(2*pi); 
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
        
        // Inner search         
        for(int jshift = -10; jshift < 11; jshift+=1) {
          for(int ishift = -10; ishift < 11; ishift+=1) {
            if(j+jshift > -1 && j+jshift < ncolx && i+ishift> -1 && i+ishift<nrowx) {
              if(land1(i+ishift,j+jshift)==1) {
                double disttoland = sqrt((jshift*jshift)+(ishift*ishift)); 
                int angletoland = 32 * (pi + atan2(jshift,ishift))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            } 
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellsz * (j + jshift)) / land2cellsz;
              int li = (land2uly - land1uly + land1cellsz * (i + ishift)) / land2cellsz;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double disttoland = sqrt((jshift*jshift)+(ishift*ishift)); 
                  int angletoland = 32 * (pi + atan2(jshift,ishift))/(2*pi); 
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
          sumfetch += minfetch(k);
          sumcos += cos(anglesector)*minfetch(k);
          sumsin += sin(anglesector)*minfetch(k);
        }
        wx32(i,j)=sumfetch;
        wx32(i,j+ncolx)=atan2(sumsin,sumcos);
        wx32(i,j+ncolx*2)=sqrt(pow((sumsin / sumfetch), 2.0) + pow((sumcos / sumfetch), 2.0));
      }
    }
    Rcout << j << std::endl;
  }
  return wx32;
}


NumericMatrix coastwx32twomap3v(IntegerMatrix land1, IntegerMatrix land2, int dx, int dwx, int cl, 
                                int nsearches, int nsectors,
                                double land1llx, double land1uly, double land1cellsz, 
                                double land2llx, double land2uly, double land2cellsz) {
  // Vector-based search. Random distance squared and random angle 
  
  int ncolx = land1.ncol();
  int nrowx = land1.nrow();
  int lncols =land2.ncol();
  int lnrows =land2.nrow();
  
  Rcout << "ncols" << ncolx << "nrows" << nrowx << std::endl;
  double pi = 3.14159265358979323846;
  
  IntegerMatrix coastcell(nrowx,ncolx);
  coastcell = isitnearcoast(land1,dx);
  
  NumericMatrix wx32(nrowx,ncolx*3);
  int vdwx = sqrt(dwx);
  
  for(int j = 0; j < ncolx; j++) {
    for(int i = 0; i < nrowx; i++) {
      if (coastcell(i,j)==cl) {
        NumericVector minfetch(nsectors);
        for(int k = 0; k < nsectors; k++) {
          minfetch(k)=dwx;
        }
        // Repeated searches along random vector - random angle and length         
        for(int m=0; m < nsearches; m++) {
          double vdist = rand() % vdwx;
          double vang = 2 * pi * ((double) rand() / RAND_MAX);
          double vx = vdist*vdist*sin(vang);
          double vy = vdist*vdist*cos(vang);
          int cx = j+vx;
          int cy = i+vy;
          //          Rcout << "search " << m << " vdist=" << vdist << " vang=" << vang << " vx=" << vx << " vy=" << vy << std::endl;
          
          if(cx > -1 && cx < ncolx && cy> -1 && cy<nrowx) {  // check if on current map
            
            if(land1(cy,cx)==1) {
              double disttoland = vdist*vdist; 
              int angletoland = nsectors * (pi + vang)/(2*pi); 
              if (angletoland> (nsectors -1)) {angletoland -=nsectors;}
              // Rcout << angletoland << "   "<< disttoland << std::endl;
              if (minfetch(angletoland)>disttoland ) {
                minfetch(angletoland)=disttoland;
              }
            }
            
          }
          else {
            // alt-map search code here
            int lj = (land1llx - land2llx + land1cellsz * (cx)) / land2cellsz;
            int li = (land2uly - land1uly + land1cellsz * (cy)) / land2cellsz;
            if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
              if(land2(li,lj)==1) {
                double disttoland = vdist*vdist; 
                int angletoland = nsectors * (pi + vang)/(2*pi); 
                if (angletoland> (nsectors -1)) {angletoland -=nsectors;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
          }
        }
        
        double sumfetch = 0;
        double sumcos = 0;
        double sumsin =0;
        for(int k = 0; k < nsectors; k++) {
          double anglesector = pi + k * 2 * pi / nsectors;
          sumfetch += minfetch(k);
          sumcos += cos(anglesector)*minfetch(k);
          sumsin += sin(anglesector)*minfetch(k);
        }
        wx32(i,j)=sumfetch;
        wx32(i,j+ncolx)=atan2(sumsin,sumcos);
        wx32(i,j+ncolx*2)=sqrt(pow((sumsin / sumfetch), 2.0) + pow((sumcos / sumfetch), 2.0));
      }  // Is it coast condition
    } // i loop
    Rcout << j << std::endl;
  } // j loop
  return wx32;
}


NumericMatrix coastwx32twomap3ve(IntegerMatrix land1, IntegerMatrix land2, int dx, int dwx, int cl, 
                                 int nsearches, int nsectors,
                                 double land1llx, double land1uly, double land1cellsz, 
                                 double land2llx, double land2uly, double land2cellsz) {
  int ncolx = land1.ncol();
  int nrowx = land1.nrow();
  int lncols =land2.ncol();
  int lnrows =land2.nrow();
  
  Rcout << "ncols " << ncolx << " nrows " << nrowx << std::endl;
  double pi = 3.14159265358979323846;
  
  IntegerMatrix coastcell(nrowx,ncolx);
  coastcell = isitnearcoast(land1,dx);
  
  NumericMatrix wx32(nrowx,ncolx*3);
  double vdwx = sqrt(sqrt(dwx));
  double dstep = vdwx * nsectors/nsearches;
  double astep = (2*pi)/nsectors;
  
  // Set up search pattern
  
  int slength = (vdwx / dstep) - 1; 
  int arrlen = slength * nsectors;
  NumericMatrix vsang(arrlen);
  NumericMatrix vdist(arrlen);
  NumericMatrix vectx(arrlen);
  NumericMatrix vecty(arrlen);
  
  Rcout << "slength " << slength << " arrlen=" << arrlen
        << " dstep " << dstep << " astep=" << astep  << std::endl;
  
  
  int counter = 1;
  for(double evend = 1; evend < vdwx; evend += dstep) {
    for(double evenang = 0; evenang < ((2*pi)-0.0001); evenang += astep) {
      vsang[counter] = evenang + astep * ((double) rand() / RAND_MAX);  // add jitter to search pattern
      vdist[counter] = pow(evend + astep * ((double) rand() / RAND_MAX), 4.0);
      vectx[counter] = 0.5 + vdist[counter]*sin(evenang);
      vecty[counter] = 0.5 + vdist[counter]*cos(evenang);
      Rcout << "search " << counter << " evend=" << evend <<
        " vdist=" << vdist[counter] << " vang=" << vsang[counter]
                  << " vx=" << vectx[counter] << " vy=" << vecty[counter] << std::endl;
        counter++;
        
    }
  }
  
  // return wx32;
  //   }
  // Cell by cell initialisation of searches
  
  for(int j = 0; j < ncolx; j++) {
    for(int i = 0; i < nrowx; i++) {
      if (coastcell(i,j)==cl) {
        NumericVector minfetch(nsectors);
        for(int k = 0; k < nsectors; k++) {
          minfetch(k)=dwx;
        }
        // Repeated searches along evenly spaced vectors - stepped angle and length
        for(int m=0; m < arrlen; m++) {
          double vang = vsang[m];
          int vx = vectx[m] ;//* (0.95 + 0.1 * ((double) rand() / RAND_MAX));
          int vy = vecty[m] ;//* (0.95 + 0.1 * ((double) rand() / RAND_MAX)); // added jitter
          int cx = j+vx;
          int cy = i+vy;
          //          Rcout << "search " << m << " vdist=" << vdist[m] << " vang=" << vang << " vx=" << vx << " vy=" << vy << std::endl;
          
          if(cx > -1 && cx < ncolx && cy> -1 && cy<nrowx) {  // check if on current map
            
            if(land1(cy,cx)==1) {
              double disttoland = vdist[m];
              int angletoland = nsectors * (pi + vang)/(2*pi);
              if (angletoland> (nsectors -1)) {angletoland -=nsectors;}
              // Rcout << angletoland << "   "<< disttoland << std::endl;
              if (minfetch(angletoland)>disttoland ) {
                minfetch(angletoland)=disttoland;
              }
            }
            
          }
          else {
            // alt-map search code here
            int lj = (land1llx - land2llx + land1cellsz * (cx)) / land2cellsz;
            int li = (land2uly - land1uly + land1cellsz * (cy)) / land2cellsz;
            if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
              if(land2(li,lj)==1) {
                double disttoland = vdist[m];
                int angletoland = nsectors * (pi + vang)/(2*pi);
                if (angletoland> (nsectors -1)) {angletoland -=nsectors;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
          }
        }  // Search pattern loop
        
        double sumfetch = 0;
        double sumcos = 0;
        double sumsin =0;
        for(int k = 0; k < nsectors; k++) {
          double anglesector = pi + k * 2 * pi / nsectors;
          sumfetch += minfetch(k);
          sumcos += cos(anglesector)*minfetch(k);
          sumsin += sin(anglesector)*minfetch(k);
        }
        wx32(i,j)=sumfetch;
        wx32(i,j+ncolx)=atan2(sumsin,sumcos);
        wx32(i,j+ncolx*2)=sqrt(pow((sumsin / sumfetch), 2.0) + pow((sumcos / sumfetch), 2.0));
      }  // Is it coast condition
    } // i loop
    Rcout << j << std::endl;
  } // j loop
  return wx32;
}

NumericMatrix coastwx32twomaplatlon3(IntegerMatrix land1, IntegerMatrix land2, int dx, int dwx, int cl, int jit,
                                     double land1llx, double land1uly, double land1cellszx, double land1cellszy, 
                                     double land2llx, double land2uly, double land2cellszx, double land2cellszy) {
  int ncolx = land1.ncol();
  int nrowx = land1.nrow();
  int lncols =land2.ncol();
  int lnrows =land2.nrow();
  
  Rcout << "ncols" << ncolx << "nrows" << nrowx << std::endl;
  double pi = 3.14159265358979323846;
  
  IntegerMatrix coastcell(nrowx,ncolx);
  coastcell = isitnearcoast(land1,dx);
  
  NumericMatrix wx32(nrowx,3*ncolx);
  
  for(int j = 0; j < ncolx; j++) {
    for(int i = 0; i < nrowx; i++) {
      if (coastcell(i,j)==cl) {
        NumericVector minfetch(32);
        for(int k = 0; k < 32; k++) {
          minfetch(k)=dwx*land1cellszy;
        }
        // Outer search        
        for(int jshift = -dwx; jshift < dwx+1; jshift+=100) {
          for(int ishift = -dwx; ishift < dwx+1; ishift+=100) {
            int jshiftj = jshift + jit*(-50 + rand() % 100);
            int ishiftj = ishift + jit*(-50 + rand() % 100);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land1(i+ishiftj,j+jshiftj)==1) {
                double yishiftj=land1cellszy*ishiftj;
                double xjshiftj=land1cellszx*jshiftj;
                double disttoland = sqrt((xjshiftj*xjshiftj)+(yishiftj*yishiftj)); 
                int angletoland = 32 * (pi + atan2(xjshiftj,yishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellszx * (j + jshiftj)) / land2cellszx;
              int li = (land2uly - land1uly + land1cellszy * (i + ishiftj)) / land2cellszy;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double yishiftj=land1cellszy*ishiftj;
                  double xjshiftj=land1cellszx*jshiftj;
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
        // Middle search
        for(int jshift = -0.1*dwx; jshift < 0.1*dwx+1; jshift+=10) {
          for(int ishift = -0.1*dwx; ishift < 0.1*dwx+1; ishift+=10) {
            int jshiftj = jshift + jit*(-5 + rand() % 10);
            int ishiftj = ishift + jit*(-5 + rand() % 10);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land1(i+ishiftj,j+jshiftj)==1) {
                double yishiftj=land1cellszy*ishiftj;
                double xjshiftj=land1cellszx*jshiftj;
                double disttoland = sqrt((xjshiftj*xjshiftj)+(yishiftj*yishiftj)); 
                int angletoland = 32 * (pi + atan2(xjshiftj,yishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellszx * (j + jshiftj)) / land2cellszx;
              int li = (land2uly - land1uly + land1cellszy * (i + ishiftj)) / land2cellszy;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double yishiftj=land1cellszy*ishiftj;
                  double xjshiftj=land1cellszx*jshiftj;
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
                // Rcout << angletoland << "   "<< disttoland << std::endl;
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
        wx32(i,j+ncolx)=atan2(sumsin,sumcos);
        wx32(i,j+ncolx*2)=sqrt(pow((sumsin / sumfetch), 2.0) + pow((sumcos / sumfetch), 2.0));
      }
    }
    Rcout << j << std::endl;
  }
  return wx32;
}



NumericMatrix coastwx32twomaplatlonft(IntegerMatrix land1, IntegerMatrix land2, int dx, int dwx, int cl, int jit,
                                      double land1llx, double land1uly, double land1cellszx, double land1cellszy, 
                                      double land2llx, double land2uly, double land2cellszx, double land2cellszy) {
  int ncolx = land1.ncol();
  int nrowx = land1.nrow();
  int lncols =land2.ncol();
  int lnrows =land2.nrow();
  
  Rcout << "ncols" << ncolx << "nrows" << nrowx << std::endl;
  double pi = 3.14159265358979323846;
  
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
        // Outer search        
        for(int jshift = -dwx; jshift < dwx+1; jshift+=100) {
          for(int ishift = -dwx; ishift < dwx+1; ishift+=100) {
            int jshiftj = jshift + jit*(-50 + rand() % 100);
            int ishiftj = ishift + jit*(-50 + rand() % 100);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land1(i+ishiftj,j+jshiftj)==1) {
                double yishiftj=land1cellszy*ishiftj;
                double xjshiftj=land1cellszx*jshiftj;
                double disttoland = sqrt((xjshiftj*xjshiftj)+(yishiftj*yishiftj)); 
                int angletoland = 32 * (pi + atan2(xjshiftj,yishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                for (int ashift = -1; ashift<2; ashift +=1) {
                  int angleref = angletoland + ashift;
                  if (angleref>31) {angleref -=32;}
                  if (angleref<0) {angleref +=32;}
                  // Rcout << "angleref" << angleref << "   "<< disttoland << std::endl;
                  if (minfetch(angleref)>disttoland ) {
                    minfetch(angleref)=disttoland;
                  }
                }
              }
            }
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellszx * (j + jshiftj)) / land2cellszx;
              int li = (land2uly - land1uly + land1cellszy * (i + ishiftj)) / land2cellszy;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double yishiftj=land1cellszy*ishiftj;
                  double xjshiftj=land1cellszx*jshiftj;
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
        // Middle search
        for(int jshift = -0.1*dwx; jshift < 0.1*dwx+1; jshift+=10) {
          for(int ishift = -0.1*dwx; ishift < 0.1*dwx+1; ishift+=10) {
            int jshiftj = jshift + jit*(-5 + rand() % 10);
            int ishiftj = ishift + jit*(-5 + rand() % 10);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land1(i+ishiftj,j+jshiftj)==1) {
                double yishiftj=land1cellszy*ishiftj;
                double xjshiftj=land1cellszx*jshiftj;
                double disttoland = sqrt((xjshiftj*xjshiftj)+(yishiftj*yishiftj)); 
                int angletoland = 32 * (pi + atan2(xjshiftj,yishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                for (int ashift = -1; ashift<2; ashift +=1) {
                  int angleref = angletoland + ashift;
                  if (angleref>31) {angleref -=32;}
                  if (angleref<0) {angleref +=32;}
                  // Rcout << "angleref midloop" << angleref << "   "<< disttoland << std::endl;
                  if (minfetch(angleref)>disttoland ) {
                    minfetch(angleref)=disttoland;
                  }
                }              }
            }
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellszx * (j + jshiftj)) / land2cellszx;
              int li = (land2uly - land1uly + land1cellszy * (i + ishiftj)) / land2cellszy;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double yishiftj=land1cellszy*ishiftj;
                  double xjshiftj=land1cellszx*jshiftj;
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
                for (int ashift = -1; ashift<2; ashift +=1) {
                  int angleref = angletoland + ashift;
                  if (angleref>31) {angleref -=32;}
                  if (angleref<0) {angleref +=32;}
                  // Rcout << "angleref innerloop" << angleref << "   "<< disttoland << std::endl;
                  if (minfetch(angleref)>disttoland ) {
                    minfetch(angleref)=disttoland;
                  }
                }
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
        // Rcout << "sumfetch" << sumfetch << std::endl;
        
        //        wx32(i,j+ncolx)=atan2(sumsin,sumcos);
        //        wx32(i,j+ncolx*2)=sqrt(pow((sumsin / sumfetch), 2.0) + pow((sumcos / sumfetch), 2.0));
      }
    }
    Rcout << j << std::endl;
  }
  return wx32;
}


NumericMatrix coastwx32twomaplatlonfmax(IntegerMatrix land1, IntegerMatrix land2, int dx, int dwx, int cl, int jit,
                                        double land1llx, double land1uly, double land1cellszx, double land1cellszy, 
                                        double land2llx, double land2uly, double land2cellszx, double land2cellszy) {
  int ncolx = land1.ncol();
  int nrowx = land1.nrow();
  int lncols =land2.ncol();
  int lnrows =land2.nrow();
  
  Rcout << "ncols" << ncolx << "nrows" << nrowx << std::endl;
  double pi = 3.14159265358979323846;
  
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
        // Outer search        
        for(int jshift = -dwx; jshift < dwx+1; jshift+=100) {
          for(int ishift = -dwx; ishift < dwx+1; ishift+=100) {
            int jshiftj = jshift + jit*(-50 + rand() % 100);
            int ishiftj = ishift + jit*(-50 + rand() % 100);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land1(i+ishiftj,j+jshiftj)==1) {
                double yishiftj=land1cellszy*ishiftj;
                double xjshiftj=land1cellszx*jshiftj;
                double disttoland = sqrt((xjshiftj*xjshiftj)+(yishiftj*yishiftj)); 
                int angletoland = 32 * (pi + atan2(xjshiftj,yishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellszx * (j + jshiftj)) / land2cellszx;
              int li = (land2uly - land1uly + land1cellszy * (i + ishiftj)) / land2cellszy;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double yishiftj=land1cellszy*ishiftj;
                  double xjshiftj=land1cellszx*jshiftj;
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
        // Middle search
        for(int jshift = -0.1*dwx; jshift < 0.1*dwx+1; jshift+=10) {
          for(int ishift = -0.1*dwx; ishift < 0.1*dwx+1; ishift+=10) {
            int jshiftj = jshift + jit*(-5 + rand() % 10);
            int ishiftj = ishift + jit*(-5 + rand() % 10);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land1(i+ishiftj,j+jshiftj)==1) {
                double yishiftj=land1cellszy*ishiftj;
                double xjshiftj=land1cellszx*jshiftj;
                double disttoland = sqrt((xjshiftj*xjshiftj)+(yishiftj*yishiftj)); 
                int angletoland = 32 * (pi + atan2(xjshiftj,yishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellszx * (j + jshiftj)) / land2cellszx;
              int li = (land2uly - land1uly + land1cellszy * (i + ishiftj)) / land2cellszy;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double yishiftj=land1cellszy*ishiftj;
                  double xjshiftj=land1cellszx*jshiftj;
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
                // Rcout << angletoland << "   "<< disttoland << std::endl;
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
        double maxfetch =0;
        double maxsin =0;
        double maxcos =0;
        for(int k = 0; k < 32; k++) {
          double anglesector = pi + k * 2 * pi / 32;
          double cosfetch = cos(anglesector)*minfetch(k)/land1cellszx;
          double sinfetch = sin(anglesector)*minfetch(k)/land1cellszy;
          double distfetch = sqrt((cosfetch*cosfetch) + (sinfetch*sinfetch));
          sumfetch += distfetch;
          sumcos += cosfetch;
          sumsin += sinfetch;
          if(distfetch>maxfetch) {
            maxfetch = distfetch;
            maxsin = sinfetch;
            maxcos = cosfetch;
          }
          
          // sumfetch += minfetch(k);
          // sumcos += cos(anglesector)*minfetch(k);
          // sumsin += sin(anglesector)*minfetch(k);
        }
        wx32(i,j)=sqrt((maxcos*maxcos) + (maxsin*maxsin));
        //        wx32(i,j)=atan2(maxsin,maxcos);
        //        wx32(i,j+ncolx)=atan2(sumsin,sumcos);
        //        wx32(i,j+ncolx*2)=sqrt(pow((sumsin / sumfetch), 2.0) + pow((sumcos / sumfetch), 2.0));
      }
    }
    Rcout << j << std::endl;
  }
  return wx32;
}

NumericMatrix coastwx32twomaplatlonomax(IntegerMatrix land1, IntegerMatrix land2, int dx, int dwx, int cl, int jit,
                                        double land1llx, double land1uly, double land1cellszx, double land1cellszy, 
                                        double land2llx, double land2uly, double land2cellszx, double land2cellszy) {
  int ncolx = land1.ncol();
  int nrowx = land1.nrow();
  int lncols =land2.ncol();
  int lnrows =land2.nrow();
  
  Rcout << "ncols" << ncolx << "nrows" << nrowx << std::endl;
  double pi = 3.14159265358979323846;
  
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
        // Outer search        
        for(int jshift = -dwx; jshift < dwx+1; jshift+=100) {
          for(int ishift = -dwx; ishift < dwx+1; ishift+=100) {
            int jshiftj = jshift + jit*(-50 + rand() % 100);
            int ishiftj = ishift + jit*(-50 + rand() % 100);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land1(i+ishiftj,j+jshiftj)==1) {
                double yishiftj=land1cellszy*ishiftj;
                double xjshiftj=land1cellszx*jshiftj;
                double disttoland = sqrt((xjshiftj*xjshiftj)+(yishiftj*yishiftj)); 
                int angletoland = 32 * (pi + atan2(xjshiftj,yishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellszx * (j + jshiftj)) / land2cellszx;
              int li = (land2uly - land1uly + land1cellszy * (i + ishiftj)) / land2cellszy;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double yishiftj=land1cellszy*ishiftj;
                  double xjshiftj=land1cellszx*jshiftj;
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
        // Middle search
        for(int jshift = -0.1*dwx; jshift < 0.1*dwx+1; jshift+=10) {
          for(int ishift = -0.1*dwx; ishift < 0.1*dwx+1; ishift+=10) {
            int jshiftj = jshift + jit*(-5 + rand() % 10);
            int ishiftj = ishift + jit*(-5 + rand() % 10);
            if(j+jshiftj > -1 && j+jshiftj < ncolx && i+ishiftj> -1 && i+ishiftj<nrowx) {
              if(land1(i+ishiftj,j+jshiftj)==1) {
                double yishiftj=land1cellszy*ishiftj;
                double xjshiftj=land1cellszx*jshiftj;
                double disttoland = sqrt((xjshiftj*xjshiftj)+(yishiftj*yishiftj)); 
                int angletoland = 32 * (pi + atan2(xjshiftj,yishiftj))/(2*pi); 
                if (angletoland>31) {angletoland -=32;}
                // Rcout << angletoland << "   "<< disttoland << std::endl;
                if (minfetch(angletoland)>disttoland ) {
                  minfetch(angletoland)=disttoland;
                }
              }
            }
            else {
              // alt-map search code here
              int lj = (land1llx - land2llx + land1cellszx * (j + jshiftj)) / land2cellszx;
              int li = (land2uly - land1uly + land1cellszy * (i + ishiftj)) / land2cellszy;
              if (lj < lncols && lj > -1 && li < lnrows && li > -1) {
                if(land2(li,lj)==1) {
                  double yishiftj=land1cellszy*ishiftj;
                  double xjshiftj=land1cellszx*jshiftj;
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
                // Rcout << angletoland << "   "<< disttoland << std::endl;
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
        double maxfetch =0;
        double maxk =0;
//        double maxcos =0;
        for(int k = 0; k < 32; k++) {
          double anglesector = pi + k * 2 * pi / 32;
          double cosfetch = cos(anglesector)*minfetch(k)/land1cellszx;
          double sinfetch = sin(anglesector)*minfetch(k)/land1cellszy;
          double distfetch = sqrt((cosfetch*cosfetch) + (sinfetch*sinfetch));
          sumfetch += distfetch;
          sumcos += cosfetch;
          sumsin += sinfetch;
          if(distfetch>maxfetch) {
            maxfetch = distfetch;
            maxk = k;
//            maxcos = cosfetch;
          }
          
        }
        //wx32(i,j)=sqrt((maxcos*maxcos) + (maxsin*maxsin));
       // wx32(i,j)=atan2(maxsin,maxcos);
        wx32(i,j)=maxk;
        //        wx32(i,j+ncolx)=atan2(sumsin,sumcos);
        //        wx32(i,j+ncolx*2)=sqrt(pow((sumsin / sumfetch), 2.0) + pow((sumcos / sumfetch), 2.0));
      }
    }
    Rcout << j << std::endl;
  }
  return wx32;
}

