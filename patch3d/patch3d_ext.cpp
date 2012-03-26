//------------------------------------------------------------
// 3d Landscape analyzer
// Calculates 3d landscape metrics based on a patch map and an elevation model
//Surface areas and distances are calculated according to:
//   Calculating Landscape Surface Area from Digital Elevation Models
//   Author(s): Jeff S. Jenness
//   Source: Wildlife Society Bulletin, Vol. 32, No. 3 (Autumn, 2004), pp. 829-839
//   Published by: Allen Press
//   Stable URL: http://www.jstor.org/stable/3784807
// Based on ideas included in the R packages:
//                     sp.surfaceArea (Barry Rowlingson)
//                     SDMTools.PatchStat(Jeremy VanDerWal)
// written by Paul Magdon 2012
//-----------------------------------------------------------
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <algorithm>
using namespace std;
//------------------------------------------------------------
typedef vector <double> data_t;
vector<data_t>::iterator p, p_end;

  struct outStructure
  {
    vector <int>    ID;
    vector <double> sarea;
    vector <double> sperimeter;
    vector <double> area;
    vector <double> perimeter;
    vector <int>    ncell;
    vector <int>    nedge;
    vector <double> para;
    vector <double> para3d;
    vector <double> frac;
    vector <double> frac3d;
    vector <double> shape;
    vector <double> shape3d;
  };


/*define the input structure*/
  struct inStructure
  {
    int ncol;
    int nrow;
    double cellsize;;
    data_t values;
  };

//-----------------------------------------------------------------
double height(double heights[], int nx, int x, int y){
  double tmp=heights[(x)+(nx)*(y)];
  return(tmp);
}
//-----------------------------------------------------------------
  void zeroit(double **array, int nrows, int ncolumns){
    int i, j;
    for(i = 0; i < nrows; i++)
      {
	for(j = 0; j < ncolumns; j++)
	  array[i][j] = 0;
      }
  }
//-----------------------------------------------------------------
inStructure extent(inStructure input){
  // create a two dimensional array tmp[Colums,Rows]
  inStructure output;
  int nx=input.ncol;
  int ny=input.ncol;
  data_t matrix=input.values;
  int size_x=nx+2;
  int size_y=ny+2;
  double **tmp;
  tmp= (double **) malloc(size_y*sizeof(double *));
  for (int i=0;i<size_y;i++){
    tmp[i]=(double *)malloc(size_x*sizeof(double));
  }

  zeroit(tmp,size_y,size_x);

  double firstline[nx];
  double lastline[nx];
  double firstcol[ny+2];
  double lastcol[ny+2];

  for(int x=0;x<nx;x++){
    for(int y=0;y<ny;y++){
      tmp[x+1][y+1]=matrix[x+(nx)*(y)];
    }
  }
     //Extract first and last line of the original data set
  for(int x=0; x<nx;x++){
    firstline[x]=tmp[x+1][1];
    lastline[x]=tmp[x+1][ny]; //Problem
  }

  //Rbind first and last line
   for(int x=0;x<nx;x++){
     tmp[x+1][0]=firstline[x];  // attach first line
     tmp[x+1][ny+1]=lastline[x]; //attach last line
   }


   //Extract first and last cols
   for(int y=0; y<ny+2;y++){
      firstcol[y]=tmp[1][y];
      lastcol[y]=tmp[nx][y];
    }

   //Cbind first and last colum
   for(int y=0;y<ny+2;y++){
     tmp[0][y]=firstcol[y];
     tmp[nx+1][y]=lastcol[y];
   }
   // Convert 2d array back to 1d array
   data_t ret;
   for(int y=0;y<(ny+2);y++){
     for(int x=0;x<(nx+2);x++){
       ret.push_back(tmp[x][y]);
     }
   }

   output.values=ret;
   output.ncol=nx+2;
   output.nrow=ny+2;
   for(int i=0;i<size_y;i++){
     free(tmp[i]);
	  }
       free(tmp);

   return(output);
}
//---------------------------------------------------------------------------------
inStructure import_ascii_grid(string path)
{

  // Declare a new vector of type data_t.
  inStructure in_data;
  data_t data=in_data.values;
  in_data.values.clear();

    ifstream infile(path.c_str());
    if(infile){
      //Importing esri asc header

      string header[6];
      for(int i=0;i<6;i++){
	getline(infile,header[i]);
      }

      string ncols,nrows,cellsizes;
      int ncol,nrow,cellsize;
      ncol=nrow=cellsize=0;
      stringstream ncolss (header[0]);
      stringstream nrowss (header[1]);
      stringstream cellsizess (header[4]);

      while (getline( ncolss, ncols, ' ' )){
	stringstream fs(ncols);
	fs >> ncol;
      }

      while (getline( nrowss, nrows, ' ' )){
	stringstream fs(nrows);
	fs >> nrow;
      }

      while (getline( cellsizess,cellsizes, ' ' )){
	stringstream fs(cellsizes);
	fs >> cellsize;
      }
      in_data.nrow=nrow;
      in_data.ncol=ncol;
      in_data.cellsize=cellsize;
      cout<<"Ncols: "<<in_data.ncol<<endl;
      cout<<"Nrows: "<<in_data.nrow<<endl;
      cout<<"Cellsize: "<<in_data.cellsize<<endl;
      while (!infile.eof()){

	string line;

	getline( infile, line );
	// now we'll use a string stream to separate the fields out of the line
	stringstream ss( line );
	string field;
	while (getline( ss, field, ' ' ))
	  {
	    //  each field is converted it to a double
	    stringstream fs( field );
	    double f = 0.0;  // (default value is 0.0)
	    fs >> f;
	    // add the newly-converted field to the end of the record
	    data.push_back( f );
	  }
      }
    }
    else
      {
	cout<<"File not found";
      }
    // Complain if something went wrong.
    if (!infile.eof())
      {
	cout << "Errors while importing raster!\n";
	//return 1;
      }

    infile.close();

    // Otherwise, list some basic information about the file.
    cout << "Raster contains " << data.size() << " records.\n";
    in_data.values=data;
    return in_data;
}
//---------------------------------------------------------------------------
int export_results(string outfile, string fileid, outStructure input){
  FILE * pFile;
  pFile =fopen (outfile.c_str(),"w");

  fprintf (pFile,"%-2s\t %-11s\t %-11s\t %-11s\t %-11s\t %-11s\t %-11s\t %-11s\t %-11s\t %-11s\t %-11s\t %-11s\t %-11s\n",
	   "ID","area3d","area","perimeter","perimeter3d","para","para3d","frac","frac3d","shape","shape3d","ncell", "nedge");
   for (int p=0;p<input.ID.size();p++){
     fprintf (pFile,"%-2d\t %-11.2f\t %-11.2f\t %-11.2f\t %-11.2f\t %-11.5f\t %-11.5f\t %-11.5f\t %-11.5f\t %-11.5f\t %-11.5f\t %-11d\t %-11d\n",
   	    input.ID[p],input.sarea[p],input.area[p],input.perimeter[p],
	      input.sperimeter[p],input.para[p],input.para3d[p],input.frac[p],input.frac3d[p],input.shape[p],input.shape3d[p],input.ncell[p],input.nedge[p]);
   }

  fclose(pFile);
  return 0;
}


//----------------------------------------------------------------------------

/* Define function for calculation of triangular areas*/
double triarea(double a, double b, double c){
  		double s;
  		s=(a+b+c)/2.0;
  		return(sqrt(s*(s-a)*(s-b)*(s-c)));
		}
//----------------------------------------------------------------------------

/* Define function for shape index calculation according to Hoechstetter et al 2008*/

double shapeindex(double perimeter,double area){
  double ret=(0.25*perimeter)/(sqrt(area));
  return ret;
}


/* Define function for the calculation of fractal dimension*/

double fractaldim(double perimeter, double area){
  double ret=(2*log(0.25*perimeter))/log(area);
  return ret;
}


//----------------------------------------------------------------------------

/*Define the patch 3d function*/
outStructure patch3d(inStructure in_ls, inStructure in_dem)
{
  int cellx=in_ls.cellsize;
  int celly=in_ls.cellsize;
  int nrow=in_ls.nrow;
  int ncol=in_ls.ncol;
  data_t ls=in_ls.values;
  outStructure output;
  int size=in_dem.ncol*in_dem.nrow;
  double *dem;
  dem= (double *) malloc(size*sizeof(double *));

  for(int i=0;i<(in_dem.nrow*in_dem.ncol);i++){
    dem[i]=in_dem.values[i];
  }

  //vector<double> IDs;
   vector<int>::iterator p, p_end;

  for(int i = 0; i<ls.size(); i++) {
    output.ID.push_back(ls[i]);
  }
  sort(output.ID.begin(),output.ID.end());             //sort vector
  p_end = unique(output.ID.begin(),output.ID.end());  // remove duplicates
  output.ID.erase(p_end,output.ID.end());

  // Set everything to zero
  for(int i=0; i<output.ID.size() ;i++){
    output.sarea.push_back(0);
    output.sperimeter.push_back(0);
    output.area.push_back(0);
    output.perimeter.push_back(0);
    output.para.push_back(0);
    output.para3d.push_back(0);
    output.frac.push_back(0);
    output.frac3d.push_back(0);
    output.shape.push_back(0);
    output.shape3d.push_back(0);
    output.ncell.push_back(0);
    output.nedge.push_back(0);
   }


  cout<< "Number of patches found: \n";
  cout <<output.ID.size()<<endl;
  int npatch = output.ID.size();
  int ncold=ncol+2;
  int nrowd=nrow+2;
  double car,cpi,np,ni;               //temporary values representing cellArea, cellPerim
  int tval, rook[4];    //values of the 9 cells of interest
  double z1,z2,z3;              // point values for each triangle
  double l1,l2,l3;              // side lengths for each triangle
  double s2 = sqrt(pow(cellx,2)+pow(celly,2)); //diagonal length depending  the raster size
  double side[]={s2,cellx,s2,celly,s2,cellx,s2,celly,s2}; // triangle side lengths
  //Rotation matrix

  int dxv[]={-1, 0, 1, 1, 1, 0,-1,-1,-1};
  int dyv[]={-1,-1,-1, 0, 1, 1, 1,0,-1};
  double l3v[]={cellx,cellx,celly,celly,cellx,cellx,celly,celly};

  cout<< "Processing col:\n"<<ncol <<endl;
  cout<< "Processing rows:\n"<<nrow<<endl;

  for (int row=0; row<nrow;row++){
    for (int col=0;col<ncol;col++){
      tval = ls[row*ncol+col];
      ni=np=0;
      car=cpi=0.0; // Set cell area and perimeter to zero
      //go clockwise through the eight neighbouring cells and collect patch ids//
      rook[0]  = (row>0)       ? ls[(row-1)*ncol+col]:-9999; // One step North
      rook[1]  = (col<ncol-1)  ? ls[(row)*ncol+(col+1)]:-9999; // One step East
      rook[2]  = (row<nrow-1)  ? ls[(row+1)*ncol+(col)]:-9999; // One step South
      rook[3]  = (col>0)       ? ls[(row)*ncol+(col-1)]:-9999; // One step West

     //cycle clockwise through the 8 trinagles and get the surface area and perimeter values
      z1 = height(dem,ncold,col+1,row+1);
      for(int tri=0;tri<8;tri++){

	// Do surface area calculations
        z2=height(dem,ncold,col+1+dxv[tri],row+1+dyv[tri]);
	z3=height(dem,ncold,col+1+dxv[tri+1],row+1+dyv[tri+1]);
        // Lenght of trinagles
	l1 = 0.5 * sqrt(side[tri]*  side[tri]  +(z1-z2)*(z1-z2));
	l2 = 0.5 * sqrt(side[tri+1]*side[tri+1]+(z1-z3)*(z1-z3));
	l3 = 0.5 * sqrt(l3v[tri]*l3v[tri]+(z2-z3)*(z2-z3));
        //calculate area of triangle
	car += triarea(l1,l2,l3);

	// Do perimeter calculations
        switch(tri){
	case 0:
	if (tval!=rook[0]){cpi+=l3;};
	break;
	case 1:
	if (tval!=rook[0]){cpi+=l3;};
	 break;
	case 2:
	 if (tval!=rook[1]){cpi+=l3;};
	   break;
	case 3:
	  if (tval!=rook[1]){cpi+=l3;};
	   break;
	case 4:
	  if (tval!=rook[2]){cpi+=l3;};
	   break;
	case 5:
	  if (tval!=rook[2]){cpi+=l3;};
	  break;
	case 6:
	  if (tval!=rook[3]){cpi+=l3;};
	  break;
	case 7:
	  if (tval!=rook[3]){cpi+=l3;};
	  break;
	}

      }
      for (int ii=0;ii<4;ii++){if (tval==rook[ii]){ni++;} else {np++;}}

      //Set cell values in output structure
      for(int ii=0;ii<npatch;ii++){
	if(output.ID[ii]==tval){
	  output.ncell[ii]      ++;
	  output.nedge[ii]      +=np;
	  output.sarea[ii]      += car;
	  output.sperimeter[ii] += cpi;
	  break;
	}
      }
    }
  }
  //Calculate 2D perimeter and area based on cell size
  for(int ii=0;ii<npatch;ii++){
    output.area[ii]=output.ncell[ii]*(cellx*celly);
    output.perimeter[ii]=output.nedge[ii]*(cellx);
    output.para[ii]=output.perimeter[ii]/output.area[ii];
    output.para3d[ii]=output.sperimeter[ii]/output.sarea[ii];
    output.frac[ii]=fractaldim(output.perimeter[ii], output.area[ii]);
    output.frac3d[ii]=fractaldim(output.sperimeter[ii],output.sarea[ii]);
    output.shape[ii]=shapeindex(output.perimeter[ii],output.area[ii]);
    output.shape3d[ii]=shapeindex(output.sperimeter[ii],output.sarea[ii]);
  }
  // Display result table
 printf("%-2s\t %-11s\t %-11s\t %-11s\t %-11s\t %-11s\t %-11s\t %-11s\t %-11s\t %-11s\t %-11s\t %-11s\t %-11s\t \n",
	"ID","Area3d","Area","Perimeter","Perimeter3d","PARA","PARA3d","FRAC","FRAC3d","SHAPE","SHAPE3d",
	   "Ncell", "Nedge");
   for (int p=0;p<npatch;p++){
     printf("%-2d\t %-11.2f\t %-11.2f\t %-11.2f\t %-11.2f\t %-11.2f\t %-11.2f\t %-11.2f\t %-11.2f\t %-11.2f\t %-11.2f\t %-11d\t %-11d\n",
   	    output.ID[p],output.sarea[p],output.area[p],output.perimeter[p],
	    output.sperimeter[p],output.para[p],output.para3d[p],output.frac[p],output.frac3d[p],output.shape[p],output.shape3d[p],output.ncell[p],output.nedge[p]);

  }
  return output;
}

int main(int argc, const char *argv[])
{
 if (argc != 4) {
   cout<<argc;
    printf("Usage: patch3d dem ls \n");
     return 1;
        }

 string dem_in=argv[1];
 string ls_in=argv[2];
 string outfile=argv[3];

 outStructure result;
 //Import elevation model
  inStructure elev=import_ascii_grid(dem_in);
  cout << "DEM Imported \n";

 //Import landscape model
 inStructure ls=import_ascii_grid(ls_in);
 cout << "Landscape Imported \n";

 // Extent elevation model by first / last row / col
 inStructure elev_ext=extent(elev);

 // Perform 3d analysis of the binary landscape
 result=patch3d(ls,elev_ext);
 // Export results
 export_results(outfile, ls_in,result);
return 0;
}
