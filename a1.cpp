#include "SImage.h"
#include "SImageIO.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include "DrawText.h"

using namespace std;

// The simple image class is called SDoublePlane, with each pixel represented as
// a double (floating point) type. This means that an SDoublePlane can represent
// values outside the range 0-255, and thus can represent squared gradient magnitudes,
// harris corner scores, etc.
//
// The SImageIO class supports reading and writing PNG files. It will read in
// a color PNG file, convert it to grayscale, and then return it to you in
// an SDoublePlane. The values in this SDoublePlane will be in the range [0,255].
//
// To write out an image, call write_png_file(). It takes three separate planes,
// one for each primary color (red, green, blue). To write a grayscale image,
// just pass the same SDoublePlane for all 3 planes. In order to get sensible
// results, the values in the SDoublePlane should be in the range [0,255].
//

// Below is a helper functions that overlays rectangles
// on an image plane for visualization purpose.

// Draws a rectangle on an image plane, using the specified gray level value and line width.
//

int FIT_MAX=0;
void overlay_rectangle(SDoublePlane &input, int _top, int _left, int _bottom, int _right, double graylevel, int width)
{
  for(int w=-width/2; w<=width/2; w++) {
    int top = _top+w, left = _left+w, right=_right+w, bottom=_bottom+w;

    // if any of the coordinates are out-of-bounds, truncate them
    top = min( max( top, 0 ), input.rows()-1);
    bottom = min( max( bottom, 0 ), input.rows()-1);
    left = min( max( left, 0 ), input.cols()-1);
    right = min( max( right, 0 ), input.cols()-1);

    // draw top and bottom lines
    for(int j=left; j<=right; j++)
	  input[top][j] = input[bottom][j] = graylevel;
    // draw left and right lines
    for(int i=top; i<=bottom; i++)
	  input[i][left] = input[i][right] = graylevel;
  }
}

// DetectedSymbol class may be helpful!
//  Feel free to modify.
//
typedef enum {NOTEHEAD=0, QUARTERREST=1, EIGHTHREST=2} Type;
class DetectedSymbol {
public:
  int row, col, width, height;
  Type type;
  char pitch;
  double confidence;
};

// Function that outputs the ascii detection output file
void  write_detection_txt(const string &filename, const vector<struct DetectedSymbol> &symbols)
{
  ofstream ofs(filename.c_str());

  for(int i=0; i<symbols.size(); i++)
    {
      const DetectedSymbol &s = symbols[i];
      ofs << s.row << " " << s.col << " " << s.width << " " << s.height << " ";
      if(s.type == NOTEHEAD)
	ofs << "filled_note " << s.pitch;
      else if(s.type == EIGHTHREST)
	ofs << "eighth_rest _";
      else
	ofs << "quarter_rest _";
      ofs << " " << s.confidence << endl;
    }
}

// Function that outputs a visualization of detected symbols
void  write_detection_image(const string &filename, const vector<DetectedSymbol> &symbols, const SDoublePlane &input)
{
  SDoublePlane output_planes[3];
  for(int i=0; i<3; i++)
    output_planes[i] = input;

  for(int i=0; i<symbols.size(); i++)
    {
      const DetectedSymbol &s = symbols[i];

      overlay_rectangle(output_planes[s.type], s.row, s.col, s.row+s.height-1, s.col+s.width-1, 255, 2);
      overlay_rectangle(output_planes[(s.type+1) % 3], s.row, s.col, s.row+s.height-1, s.col+s.width-1, 0, 2);
      overlay_rectangle(output_planes[(s.type+2) % 3], s.row, s.col, s.row+s.height-1, s.col+s.width-1, 0, 2);

      if(s.type == NOTEHEAD)
	{
	  char str[] = {s.pitch, 0};
	  draw_text(output_planes[0], str, s.row, s.col+s.width+1, 0, 2);
	  draw_text(output_planes[1], str, s.row, s.col+s.width+1, 0, 2);
	  draw_text(output_planes[2], str, s.row, s.col+s.width+1, 0, 2);
	}
    }

  SImageIO::write_png_file(filename.c_str(), output_planes[0], output_planes[1], output_planes[2]);
}

void one_d_convolve(const SDoublePlane &matrixone,const SDoublePlane &matrixtwo)
{

	  double sum=0;
	  int counter=0;

	  if(matrixone.cols()>matrixtwo.cols())
		  counter=matrixone.cols();
	  else
		  counter=matrixtwo.cols();
	  for(int i=0;i<counter;i++)
	      {
		  if(i<matrixone.cols() and i<matrixtwo.cols())
		  sum=sum+matrixone[0][i]*matrixtwo[0][i];
	      }
	  cout<<sum;
}

SDoublePlane one_d_row_convolve(const SDoublePlane &matrixone,const SDoublePlane &matrixtwo)
{
	SDoublePlane res_row_convolve(1,matrixone.cols());
	for(int i=0;i<matrixone.cols();i++)
	{
		for(int j=0;j<matrixone.rows();j++)
		{
			if(j<matrixtwo.rows())
			{
			res_row_convolve[0][i]=res_row_convolve[0][i]+matrixone[j][i]*matrixtwo[j][0];
			}
		}
	}
	return res_row_convolve;
}

SDoublePlane one_d_col_convolve(const SDoublePlane &matrixone,const SDoublePlane &matrixtwo)
{
	SDoublePlane res_col_convolve(matrixone.rows(),1);
	for(int i=0;i<matrixone.rows();i++)
	{
		for(int j=0;j<matrixone.cols();j++)
		{
			if(j<matrixtwo.cols())
			{
			res_col_convolve[0][i]=res_col_convolve[0][i]+matrixone[i][j]*matrixtwo[0][j];
			}
		}
	}
	return res_col_convolve;
}



// The rest of these functions are incomplete. These are just suggestions to
// get you started -- feel free to add extra functions, change function
// parameters, etc.

// Convolve an image with a separable convolution kernel
//
SDoublePlane convolve_separable(const SDoublePlane &input, const SDoublePlane &row_filter, const SDoublePlane &col_filter)
{
  SDoublePlane output(input.rows(), input.cols());
	int filter_center_row=row_filter.rows()/2;
	  int filter_center_col=col_filter.cols()/2;
	  SDoublePlane intermediated_array(1,input.cols());
	       SDoublePlane sub_array(input.rows(),input.cols());

      for(int i=0;i<input.rows();i++)
      {
    	  for(int j=0;j<input.cols();j++)
    	  {
    		  int a=0;
    		  for(int k=i-filter_center_row;k<i-filter_center_row+row_filter.rows();k++)
    		  {
    			  int b=0;
    			  for(int l=j-filter_center_col;l<j-filter_center_col+col_filter.cols();l++)
    			  {
    				  if(k>=0 and l>=0 and k<input.rows() and l<input.cols())
    				  {
    					  sub_array[a][b]=input[k][l];
    				  }
    				  b=b+1;
    			  }
    			  a=a+1;
    		  }
    		  intermediated_array=one_d_row_convolve(sub_array,row_filter);
    		  output[i][j]=one_d_col_convolve(intermediated_array,col_filter)[0][0];
    	  }
      }
  return output;
}

// Convolve an image with a separable convolution kernel
//
SDoublePlane convolve_general(const SDoublePlane &input, const SDoublePlane &filter)
{
  SDoublePlane output(input.rows(), input.cols());
  int filter_center_row=filter.rows()/2;
  int filter_center_col=filter.cols()/2;

  for(int i=0;i<input.rows();i++)
  {
	  for(int j=0;j<input.cols();j++)
	  {
		  int sum=0;
		  for(int k=0;k<filter.rows();k++)
		  {
			  for(int l=0;l<filter.cols();l++)
			  {
				  if(i-filter_center_row+k>=0 and i-filter_center_row+k<input.rows() and j-filter_center_col+l>=0 and j-filter_center_col+l<input.cols())
				  {
						  sum=sum+input[i-filter_center_row+k][j-filter_center_col+l]*filter[k][l];
				  }
			  }
		  }
		  output[i][j]=sum;
	  }
  }


  return output;
}

SDoublePlane gaussianBlurkernel(float sigma)
{
  SDoublePlane kernel(21,21);
  int ac, bc;
  float g, sum;
  sum = 0;
  for(int a=0; a<kernel.rows(); a++) {
    for(int b=0; b<kernel.cols(); b++) {

      ac = a - (kernel.rows()-1)/2;
      bc = b - (kernel.cols()-1)/2;

      g = exp(-(ac*ac+bc*bc)/(2*sigma*sigma));
      sum += g;
      kernel[a][b] = g;
    }
  }
  // Normalization
    for(int a=0; a<kernel.rows(); a++) {
      for(int b=0; b<kernel.cols(); b++) {
        kernel[a][b] /= sum;
      }
    }
    return kernel;
}



// Apply a sobel operator to an image, returns the result
//
SDoublePlane sobel_gradient_filter(const SDoublePlane &input, bool _gx)
{
  SDoublePlane output(input.rows(), input.cols());

  // Implement a sobel gradient estimation filter with 1-d filters


  return output;
}

SDoublePlane sobel_edge_detector(const SDoublePlane &input,const SDoublePlane &xfilter,const SDoublePlane &yfilter )
{
  SDoublePlane output(input.rows(), input.cols());
  SDoublePlane x_output(input.rows(), input.cols());
  SDoublePlane y_output(input.rows(), input.cols());
  x_output=convolve_general(input,xfilter);
  y_output=convolve_general(input,yfilter);
  for(int i=0;i<output.rows();i++)
  {
	  for(int j=0;j<output.cols();j++)
	  {
		  output[i][j]=sqrt(pow (x_output[i][j],2)+pow (y_output[i][j],2)) ;
		  //output[i][j]=abs(x_output[i][j])+abs (y_output[i][j]) ; //Could be used for faster computation
	  }
  }

  // Implement an edge detector of your choice, e.g.
  // use your sobel gradient operator to compute the gradient magnitude and threshold

  return output;
}

// Apply an edge detector to an image, returns the binary edge map
//
SDoublePlane find_edges(const SDoublePlane &input, double thresh)
{
  SDoublePlane output(input.rows(), input.cols());
  for(int i=0;i<input.rows();i++)
    {
  	  for(int j=0;j<input.cols();j++)
  	  {
  		if(input[i][j]>thresh)
  			input[i][j]=255;
  		     // else if(input[i][j]<(thresh/2))
  			//input[i][j]=0;
  		else
  		    	input[i][j]=0;
  		output[i][j]=255-(unsigned char)(input[i][j]);
  	  }
    }


  // Implement an edge detector of your choice, e.g.
  // use your sobel gradient operator to compute the gradient magnitude and threshold

  return output;
}

SDoublePlane output_edge_detector(const SDoublePlane &input)
{
	SDoublePlane output(input.rows(), input.cols());
	SDoublePlane gaussian_kernel(21,21);
	  	SDoublePlane sobel_x_filter(3,3);
	  	SDoublePlane sobel_y_filter(3,3);
	  	SDoublePlane gaussianBlurImage(input.rows(),input.cols());
	  	float sigma=0.05;// Sigma can be - 0.01, 1.0, 5.0, and 10.0
	  	double thresh;//Could be set to 200, if not dynamic
	  	int x_filter[3][3]={
	  	    		  {-1,0,+1,},
	  	    		  {-2,0,2,},
	  	    		  {-1,0,1}
	  	      };
	  	  int y_filter[3][3]={
	  			  	  {1,2,+1,},
	  			  	  {0,0,0,},
	  			  	  {-1,-2,-1}
	  	        };
	  	for(int i=0; i<sobel_x_filter.rows(); i++)
	  	  {
	  		  for(int j=0; j<sobel_x_filter.cols(); j++)
	  		    {
	  			sobel_x_filter[i][j] = x_filter[i][j];
	  		    }
	  	  }
	  	for(int i=0; i<sobel_y_filter.rows(); i++)
	  	  {
	  		  for(int j=0; j<sobel_y_filter.cols(); j++)
	  		    {
	  			sobel_y_filter[i][j] = y_filter[i][j];
	  		    }
	  	  }
	  	gaussian_kernel=gaussianBlurkernel(sigma);
	  	gaussianBlurImage=convolve_general(input,gaussian_kernel);
	  	output=sobel_edge_detector(input,sobel_x_filter,sobel_y_filter);
	  	for (int i = 0; i < output.rows(); i++)
	  	    {
	  		for (int j = 0; j < output.cols(); j++)
	  		  	    {
	  	        if (output[i][j] >= thresh)
	  	        {
	  	        	thresh = output[i][j];
	  	        }
	  	    }
	  	    }


	  	output=find_edges(output,thresh/5);
	  	return output;

}

SDoublePlane binary_image(const SDoublePlane &input)
{
	SDoublePlane output(input.rows(), input.cols());

	for(int i=0;i<input.rows();i++)
	    {
	  	  for(int j=0;j<input.cols();j++)
	  	  {
	  		  if(input[i][j]==255)
	  		  {
	  			 output[i][j]=1;
	  		  }
	  		  else
	  		  {
	  			output[i][j]= 0;
	  		  }
	  	  }
	    }

	return output;
}

SDoublePlane closet_edge_pixel(const SDoublePlane &input_binary_image)
{
	SDoublePlane output(input_binary_image.rows(), input_binary_image.cols());
	SDoublePlane gamma(input_binary_image.rows(), input_binary_image.cols());
	for(int i=0;i<input_binary_image.rows();i++)
			    {
			  	  for(int j=0;j<input_binary_image.cols();j++)
			  	  {
			  		if(input_binary_image[i][j]==0)
			  			  		  {
			  						gamma[i][j]= 999999999;
			  			  		  }
			  			  		  else
			  			  		  {
			  			  			gamma[i][j]= 0;
			  			  		  }
			  	  }
			    }



	for(int i=0;i<input_binary_image.rows();i++)
		    {
		  	  for(int j=0;j<input_binary_image.cols();j++)
		  	  {
		  		vector<float> row(input_binary_image.rows());
		  		  for(int k=0;k<input_binary_image.rows();k++)
		  		  {

		  			vector<float> colm(input_binary_image.cols());
		  			  for(int l=0;l<input_binary_image.cols();l++)
		  				  	  {
		  				  	  	  if(gamma[k][l]==0)
		  				  	  	  {
		  				  	  	  if(k<gamma.rows() and l<gamma.cols())
		  				  	  	  {
		  				  		  float tempuu=sqrt(pow(i-k ,2)+pow(j-l ,2));
		  				  		  float dis=gamma[k][l] + tempuu;
		  				  		  colm[l]=dis;
		  				  		//output[i][j]= gamma[k][l] + sqrt(pow(i-k ,2)+pow(j-l ,2));
		  				  	  	  }
		  				  	  	  }
		  				  	  	  else
		  				  	  	  {
		  				  	  		  colm[l]=999999999;
		  				  	  	  }

		  				  	  }
		  			int mic = colm[0];
		  			    for(int i=0;i<colm.size();i++)
		  			    {
		  			        if(colm[i]<mic)
		  			        mic=colm[i];
		  			    }
		  			  row[k]=mic;
		  		  }
		  		int mic = row[0];
		  				  			    for(int i=0;i<row.size();i++)
		  				  			    {
		  				  			        if(row[i]<mic)
		  				  			        mic=row[i];
		  				  			    }
		  		output[i][j]=mic;

		  				    }
		  	  }



	return output;
}

SDoublePlane covert_image(const SDoublePlane &input_image)
{
	SDoublePlane output(input_image.rows(), input_image.cols());
	for(int i=0; i<input_image.rows();i++)
	{
		for(int j=0;j<input_image.cols();j++)
		{
			if(input_image[i][j]==0)
			{
				output[i][j]=1;
			}
			else
			{
				output[i][j]=0;
			}
		}
	}
	return output;
}

SDoublePlane FIT(const SDoublePlane &input_image,const SDoublePlane &template_image,vector<DetectedSymbol> &symbols,Type type)
{

	SDoublePlane output(input_image.rows(), input_image.cols());
	SDoublePlane closest_edge_output(input_image.rows(), input_image.cols());
	SDoublePlane input_binary_image=binary_image(input_image);

		SDoublePlane template_binary_image=binary_image(template_image);

		closest_edge_output=closet_edge_pixel(input_binary_image);




				for(int i=0;i<closest_edge_output.rows();i++)
				    {
				  	  for(int j=0;j<closest_edge_output.cols();j++)
				  	  {
				  		float tempu=0.0;
				  		for(int k=0;k<template_binary_image.rows()-1;k++)
				  		{
				  			for(int l=0;l<template_binary_image.cols()-1;l++)
				  				{
				  					if(i+k <closest_edge_output.rows() and j+l <closest_edge_output.cols())
				  					{

				  						tempu=tempu+(template_binary_image[k][l])*(closest_edge_output[i+k][j+l]);

				  					}

				  				}

				  		 }

				  		output[i][j]=tempu;
				  	  }
				    }


				int dontknow=0;
				for(int i=0;i<output.rows();i++){
										for(int j=0;j<output.cols();j++){
											if(output[i][j]>dontknow)
											{
												dontknow=output[i][j];
											}
										}
				}

				for(int i=0;i<output.rows();i++){
						for(int j=0;j<output.cols();j++){
							if(output[i][j] <dontknow/17){
								i=i+template_binary_image.rows();
								j=j+template_binary_image.cols();
									  DetectedSymbol s;
								      s.row =i-template_binary_image.rows()/3;
								      s.col =j-template_binary_image.cols()/3;
								      s.width = template_binary_image.cols()-4;
								      s.height = template_binary_image.rows()-4;
								      s.type = type;
								      symbols.push_back(s);

							}


						}

				}

	return output;
}


//
// This main file just outputs a few test images. You'll want to change it to do
//  something more interesting!
//
int main(int argc, char *argv[])
{
  if(!(argc >= 2))
    {
      cerr << "usage: " << argv[0] << " input_image" << endl;
      return 1;
    }

  string input_filename(argv[1]);
  SDoublePlane input_image= SImageIO::read_png_file(input_filename.c_str());
  // test step 2 by applying mean filters to the input image
  SDoublePlane mean_filter(3,3);
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      mean_filter[i][j] = 1/9.0;

//-----------------------------------------------------------------------------------------------------
/*
 * Convolution Functionality below by debasis dwivedy
 *
 */

  SDoublePlane output_image(input_image.rows(),input_image.cols());
  SDoublePlane gen_convolve_filter(3,3);
  SDoublePlane sep_convolve_row_filter(3,1);
  SDoublePlane sep_convolve_col_filter(1,3);
  /*
   * General Convolution
   *
   */
  int kernel[3][3]={
  		  {-1, -2, -1,},
  		  {0, 0, 0,},
  		  {1, 2, 1}
    };
  int new_kernel[3][3]={0};
    for (int i=0;i<3;i++)
  		{
  	  for(int j=0;j<3;j++)
  	  {
  		  new_kernel[2-i][2-j]=kernel[i][j];
  	  }
  		}
  for(int i=0; i<3; i++)
      for(int j=0; j<3; j++)
    	  gen_convolve_filter[i][j] = new_kernel[i][j];
  //output_image = convolve_general(input_image, gen_convolve_filter);

/*
 * Separable Convolution
 *
 */

  int row_filter[3][1]={
    		  {1,},
    		  {2,},
    		  {1}
      };
  int col_filter[1][3]={
      		  {1, 2, 1}
        };
  for(int i=0; i<sep_convolve_row_filter.rows(); i++)
  {
	  for(int j=0; j<sep_convolve_row_filter.cols(); j++)
	    {
      	  sep_convolve_row_filter[i][j] = row_filter[i][j];
	    }
  }

  for(int i=0; i<sep_convolve_col_filter.rows(); i++)
  {
  for(int j=0; j<sep_convolve_col_filter.cols(); j++)
  {
      	  sep_convolve_col_filter[i][j] = col_filter[i][j];
  }
  }


/*
 * End of Convolution
 *
 */
//----------------------------------------------------------------------------------------------------------
/*
 *
 * Sobel gradient and edge detector
 *
 * To get the right gradient value divide sobel_x_filter and sobel_y_filter by 1/8
 */

  string template1_filename(argv[2]);
      SDoublePlane template1_image= SImageIO::read_png_file(template1_filename.c_str());
      string template2_filename(argv[3]);
          SDoublePlane template2_image= SImageIO::read_png_file(template2_filename.c_str());
          string template3_filename(argv[4]);
              SDoublePlane template3_image= SImageIO::read_png_file(template3_filename.c_str());
      SDoublePlane edge_input_image(input_image.rows(),input_image.cols());
        SDoublePlane edge_template_image1(template1_image.rows(),template1_image.cols());
        SDoublePlane edge_template_image2(template2_image.rows(),template2_image.cols());
        SDoublePlane edge_template_image3(template3_image.rows(),template3_image.cols());
        edge_input_image=output_edge_detector(input_image);
        edge_template_image1=output_edge_detector(template1_image);
        edge_template_image2=output_edge_detector(template2_image);
        edge_template_image3=output_edge_detector(template3_image);


/*
 *
 * End of Sobel gradient and edge detector
 */
//----------------------------------------------------------------------------------------------------------------
  // randomly generate some detected symbols -- you'll want to replace this
  //  with your symbol detection code obviously!
  vector<DetectedSymbol> symbols;

      vector<DetectedSymbol> symbol;
      write_detection_image("edge_template_image1.png", symbol, edge_template_image1);
      write_detection_image("edge_template_image2.png", symbol, edge_template_image2);
      write_detection_image("edge_template_image3.png", symbol, edge_template_image3);
      write_detection_image("edges.png", symbol, edge_input_image);

//Code for cross correlation

  FIT(edge_input_image,edge_template_image1,symbols,NOTEHEAD);

 FIT(edge_input_image,template2_image,symbols,QUARTERREST);
  FIT(edge_input_image,edge_template_image3,symbols,EIGHTHREST);
  write_detection_image("detected5.png", symbols, input_image);
}
