#include <SImage.h>
#include <SImageIO.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <DrawText.h>

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



// The rest of these functions are incomplete. These are just suggestions to 
// get you started -- feel free to add extra functions, change function
// parameters, etc.

// Convolve an image with a separable convolution kernel
//
SDoublePlane convolve_separable(const SDoublePlane &input, const SDoublePlane &row_filter, const SDoublePlane &col_filter)
{
  SDoublePlane output(input.rows(), input.cols());
  SDoublePlane output_temp(input.rows(), input.cols());


  // Convolution code here
  	int k=col_filter.rows();
    int c=col_filter.rows()/2;

    int imageRows=input.rows();
    int imageCols=input.cols();

    //using row filter
    for(int i=c;i<imageRows-c;i++){
  	  for(int j=c;j<imageCols-c;j++){
  		  double temp=0;

  		  	  for(int v=0;v<k;v++){
  		  			temp=temp+row_filter[0][k-1-v] * input[i-c][j+k-1-v-c];

  		  		  }

  			output_temp[i][j]=temp;
  	  }

    }
    //using column filter
    for(int i=c;i<imageRows-c;i++){
      	  for(int j=c;j<imageCols-c;j++){
      		  double temp=0;

      		  	  for(int u=0;u<k;u++){
      		      		 temp=temp+col_filter[k-1-u][0] * output_temp[i+k-1-u-c][j-c];

      		      		  }

      			output[i][j]=temp;
      	  }

        }


    return output;
}

// Convolve an image with a separable convolution kernel


SDoublePlane convolve_general(const SDoublePlane input, const SDoublePlane filter)
{
  SDoublePlane output(input.rows(), input.cols());

  // Convolution code here
  int k=filter.rows();
  int c=k/2;

  int imageRows=input.rows();
  int imageCols=input.cols();

  int m,n;
  for(int i=c;i<imageRows-c;i++){
	  for(int j=c;j<imageCols-c;j++){
		  double temp=0;
		  for(int u=0;u<k;u++){
			  	  	  m=k-1-u;
		  			  for(int v=0;v<k;v++){
		  				n=k-1-v;
		  				temp=temp+filter[m][n] * input[i+m-c][j+n-c];

		  			  }

		  		  }

		  output[i][j]=temp;
	  }
  
  }
  return output;
}
/*
 * This method binarize the image to be used for hamming distance computation
 */
SDoublePlane binarize_image(const SDoublePlane image){
	SDoublePlane output_image(image.rows(),image.cols());
	for(int i=0;i<image.rows();i++){
		for(int j=0;j<image.cols();j++){
			if(image[i][j]<=150){
				output_image[i][j]=0;
			}
			else{
				output_image[i][j]=1;
			}
		}
	}
	return output_image;
}
/*
 * This is a helper method to display the pixel values of an image.
 */
void display_pixel_values(const SDoublePlane image){
	int min=9999;
	for(int i=0;i<image.rows();i++){
		for(int j=0;j<image.cols();j++){
			cout<< image[i][j]<<" ";
		}
		cout<< endl;
	}
}


void block_image_area(SDoublePlane &hamming_matrix,int i,int j,int template_width,int template_height){
	for(int m=i;m<i+template_height;m++){
		for(int n=j;n<j+template_width;n++){
			if(m!=i && n!=j){
				hamming_matrix[m][n]=-1;
			}
		}
	}
}
/*
 * This method detects the symbols.
 */
void detect_symbols(const SDoublePlane input,const SDoublePlane note_template,vector<DetectedSymbol> &symbols,Type type){

	SDoublePlane hamming_matrix(input.rows(),input.cols());

// Binarize the input image and the template
	SDoublePlane binarized_input_image=binarize_image(input);
	SDoublePlane binarized_template_image=binarize_image(note_template);

// Hamming Distance Computation
	double temp;
	int max=-1;
	int c=note_template.cols()/2;
	int d=note_template.rows()/2;
	int m,n;
	for(int i=d;i<binarized_input_image.rows()-d;i++){
		for(int j=c;j<binarized_input_image.cols()-c;j++){
			temp=0;
			for(int k=0;k<binarized_template_image.rows();k++){
				m=binarized_template_image.rows()-1-k;
				for(int l=0;l<binarized_template_image.cols();l++){
					n=binarized_template_image.cols()-1-l;
					temp=temp+binarized_input_image[i+m-d][j+n-c] * binarized_template_image[m][n]+(1-binarized_input_image[i+m-d][j+n-c]) * (1-binarized_template_image[m][n]);
				}
			}
			hamming_matrix[i][j]=temp;
			if(temp>max)
				max=temp;
		}
	}
	for(int i=0;i<hamming_matrix.rows();i++){
		for(int j=0;j<hamming_matrix.cols();j++){
			if(hamming_matrix[i][j] >0 && hamming_matrix[i][j]>=0.91*max){
					  DetectedSymbol s;
				      s.row =i-binarized_template_image.rows()/2;
				      s.col =j-binarized_template_image.cols()/2;
				      s.width = binarized_template_image.cols();
				      s.height = binarized_template_image.rows();
				      s.type = type;
				      symbols.push_back(s);
				      block_image_area(hamming_matrix,s.row,s.col,s.width,s.height);
			}
		}
	}




}
/*
 * This method detect the co-ordinates of the staff lines in the image using a 3X3 kernel
 * |-1,-1,-1|
 * | 2, 2, 2|
 * |-1,-1,-1|
 */
vector<int> detect_lines(const SDoublePlane input) {
	SDoublePlane output_lines(input.rows(), input.cols());
	SDoublePlane horizontal_line_kernel(3, 3);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if (i == 0 || i == 2)
				horizontal_line_kernel[i][j] = -1;
			else
				horizontal_line_kernel[i][j] = 2;
		}
	}
		int m,n;
		for (int i = 1; i < input.rows() - 1; i++) {
			for (int j = 1; j < input.cols() - 1; j++) {
				int temp = 0;
				for (int u = 0; u < 3; u++) {
					m=2-u;
					for (int v = 0; v < 3; v++) {
						n=2-v;
						temp = temp
								+ input[i + m - 1][j + n - 1]
										* horizontal_line_kernel[m][n];
					}
				}
				if (temp<-460)
					output_lines[i][j] = 0;
				else
					output_lines[i][j] = 255;

			}
		}



	vector<int> indices;
	  for(int i=0;i<output_lines.rows();i++){
		  int count_black=0;
		  for(int j=0;j<output_lines.cols();j++){
			  if(output_lines[i][j]==0){
				  count_black++;
			  }
		  }
		  if(count_black>0.20*output_lines.cols()){

			  indices.push_back(i);

			  }


	  }
return indices;

}

/*
 * This method computes the pitch of the detected notes.
 */
void compute_pitch(vector<int> line_positions,vector<DetectedSymbol> &symbols){

for(int i=0;i<symbols.size();i++){
	for(int j=0;j<=(line_positions.size()/10-1);j++){
		if(symbols[i].type==NOTEHEAD){
			if(symbols[i].row+5 <line_positions[1+j*10] || (symbols[i].row+5>line_positions[4+j*10]-3 && symbols[i].row+5<line_positions[4+j*10]+3) ||(symbols[i].row+5>line_positions[6+j*10]+5 && symbols[i].row+5<line_positions[6+j*10]+8)||(symbols[i].row+5>line_positions[10+j*10]-3 && symbols[i].row+5<line_positions[10+j*10]+3) ){
							symbols[i].pitch='G';
							break;
						}
						else if((symbols[i].row+5>line_positions[1+j*10]-3 && symbols[i].row+5 <=line_positions[1+j*10]+3)||(symbols[i].row+5>line_positions[4+j*10]+5 && symbols[i].row+5 <=line_positions[4+j*10]+8) ||(symbols[i].row+5>line_positions[7+j*10]-3 && symbols[i].row+5 <=line_positions[7+j*10]+3)||(symbols[i].row+5>line_positions[10+j*10]+5 && symbols[i].row+5 <=line_positions[10+j*10]+8)){
							symbols[i].pitch='F';
							break;
						}
						else if((symbols[i].row+5>line_positions[1+j*10]+5&& symbols[i].row+5 <=line_positions[1+j*10]+8)||(symbols[i].row+5>line_positions[5+j*10]-3 && symbols[i].row+5 <=line_positions[5+j*10]+3)||(symbols[i].row+5>line_positions[7+j*10]+5 && symbols[i].row+5 <=line_positions[7+j*10]+8)||((symbols[i].row+5>line_positions[10+j*10]+8 && symbols[i].row+5 <=line_positions[10+j*10]+14))){
							symbols[i].pitch='E';
							break;
						}
						else if((symbols[i].row+5>line_positions[2+j*10]-3 && symbols[i].row+5 <=line_positions[2+j*10]+3)||(symbols[i].row+5>line_positions[5+j*10]+5&& symbols[i].row+5 <=line_positions[5+j*10]+8)||(symbols[i].row+5>line_positions[8+j*10]-3&& symbols[i].row+5 <=line_positions[8+j*10]+3)){
								symbols[i].pitch='D';
								break;
							}
						else if((symbols[i].row+5>line_positions[2+j*10]+5 && symbols[i].row+5 <=line_positions[2+j*10]+8)||(symbols[i].row+5>line_positions[5+j*10]+8 && symbols[i].row+5 <=line_positions[5+j*10]+10)||(symbols[i].row+5>line_positions[6+j*10]-14 && symbols[i].row+5 <=line_positions[6+j*10]-8)||(symbols[i].row+5>line_positions[8+j*10]+5 && symbols[i].row+5 <=line_positions[8+j*10]+8)){
									symbols[i].pitch='C';
									break;
								}
						else if((symbols[i].row+5>line_positions[3+j*10]-3 && symbols[i].row+5 <=line_positions[3+j*10]+3)||(symbols[i].row+5>line_positions[5+j*10]+16 && symbols[i].row+5 <line_positions[5+j*10]+20)||(symbols[i].row+5>line_positions[6+j*10]-9 && symbols[i].row+5 <line_positions[6+j*10]-3)||(symbols[i].row+5>line_positions[9+j*10]-3 && symbols[i].row+5 <line_positions[9+j*10]+3)){
										symbols[i].pitch='B';
										break;
									}
						else if((symbols[i].row+5>line_positions[3+j*10]+5 && symbols[i].row+5 <=line_positions[3+j*10]+8) || (symbols[i].row+5>line_positions[6+j*10]-3 && symbols[i].row+5 <=line_positions[6+j*10]+3)||(symbols[i].row+5>line_positions[9+j*10]+5 && symbols[i].row+5 <=line_positions[9+j*10]+8)){
											symbols[i].pitch='A';
											break;
										}

			}

	}
}
}


// Apply a sobel operator to an image, returns the result
// 
SDoublePlane sobel_gradient_filter(const SDoublePlane &input, bool _gx)
{
  SDoublePlane output(input.rows(), input.cols());

  // Implement a sobel gradient estimation filter with 1-d filters
  

  return output;
}

// Apply an edge detector to an image, returns the binary edge map
// 
SDoublePlane find_edges(const SDoublePlane &input, double thresh=0)
{
  SDoublePlane output(input.rows(), input.cols());

  // Implement an edge detector of your choice, e.g.
  // use your sobel gradient operator to compute the gradient magnitude and threshold
  
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
  int filtersize=3;
  SDoublePlane mean_filter(filtersize,filtersize);
  for(int i=0; i<filtersize; i++)
    for(int j=0; j<filtersize; j++)
      mean_filter[i][j] = 1/((double) filtersize*filtersize);
  SDoublePlane output_1 = convolve_general(input_image, mean_filter);
//
  // Separable filter
  SDoublePlane row_filter(1,filtersize);
  for(int i=0; i<filtersize; i++)
	  row_filter[0][i] = 1/((double) filtersize);
  SDoublePlane col_filter(filtersize,1);
  for(int j=0; j<filtersize; j++)
  	  col_filter[j][0] = 1/((double) filtersize);
  SDoublePlane output_2 = convolve_separable(input_image, row_filter,col_filter);

  string template1_filename(argv[2]);
  SDoublePlane template1_image= SImageIO::read_png_file(template1_filename.c_str());

  string template2_filename(argv[3]);
  SDoublePlane template2_image= SImageIO::read_png_file(template2_filename.c_str());

  string template3_filename(argv[4]);
  SDoublePlane template3_image= SImageIO::read_png_file(template3_filename.c_str());

//
//  // randomly generate some detected symbols -- you'll want to replace this
//  //  with your symbol detection code obviously!
  vector<DetectedSymbol> symbols;

  detect_symbols(output_2,template1_image,symbols,NOTEHEAD);
  detect_symbols(output_2,template2_image,symbols,QUARTERREST);
  detect_symbols(output_2,template3_image,symbols,EIGHTHREST);

  vector<int> line_positions=detect_lines(input_image);
  compute_pitch(line_positions,symbols);





  write_detection_txt("detected.txt", symbols);
  write_detection_image("scores4_convolve_genereal.png", symbols, output_1);
  write_detection_image("scores4.png", symbols, output_2);
  write_detection_image("detected4.png", symbols, input_image);


}
