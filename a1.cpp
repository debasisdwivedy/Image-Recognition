#include "SImage.h"
#include "SImageIO.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include "DrawText.h"
#include "opencv2/opencv.hpp"

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

  for(int i=0; i<symbols.size(); i++){
    const DetectedSymbol &s = symbols[i];

    overlay_rectangle(output_planes[s.type], s.row, s.col, s.row+s.height-1, s.col+s.width-1, 255, 2);
    overlay_rectangle(output_planes[(s.type+1) % 3], s.row, s.col, s.row+s.height-1, s.col+s.width-1, 0, 2);
    overlay_rectangle(output_planes[(s.type+2) % 3], s.row, s.col, s.row+s.height-1, s.col+s.width-1, 0, 2);

    if(s.type == NOTEHEAD){
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

  int filter_center_row=filter.rows()/2;
  int filter_center_col=filter.cols()/2;

  for(int i=0;i<input.rows();i++){
    for(int j=0;j<input.cols();j++){
      int sum=0;
      for(int k=0;k<filter.rows();k++){
        for(int l=0;l<filter.cols();l++){
          if(i-filter_center_row+k>=0 and i-filter_center_row+k<input.rows() and j-filter_center_col+l>=0 and j-filter_center_col+l<input.cols()){
            sum=sum+input[i-filter_center_row+k][j-filter_center_col+l]*filter[k][l];
          }
        }
      }
      output[i][j]=sum;
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


/**
 * This method prevents multiple overlapping squares to be drawn around a symbol
 */
void block_image_area(SDoublePlane &hamming_matrix,int symbol_x_coordinate,int symbbol_y_coordinate,int template_width,int template_height){
	for(int m=symbol_x_coordinate;m<symbol_x_coordinate+template_height;m++){
		for(int n=symbbol_y_coordinate;n<symbbol_y_coordinate+template_width;n++){
			if(m!=symbol_x_coordinate && n!=symbbol_y_coordinate){
				hamming_matrix[m][n]=-1;
			}
		}
	}
}
/*
 * This method detects the symbols by using Hamming Distance as the similarity function.
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
* For Q4, this method detect the coordinates of the staff lines in the image by convolution of the below kernel with the input image.
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
 * This method computes the pitch of the detected note heads.
 * Params:
 * line_positions: row coordinates of detected staff lines
 * symbols: coordinates of NOTEHEADS detected
 * avg_scale: height of a NOTEHEAD
 */
void compute_pitch(vector<int> line_positions,vector<DetectedSymbol> &symbols,int avg_scale){

for(int i=0;i<symbols.size();i++){

	for(int j=0;j<=(line_positions.size()/10-1);j++){

		if(symbols[i].type==NOTEHEAD){
					if(symbols[i].row+avg_scale/2 <line_positions[0+j*10] || (symbols[i].row+avg_scale/2>line_positions[3+j*10]-avg_scale/3 && symbols[i].row+avg_scale/2<line_positions[3+j*10]+avg_scale/3) ||(symbols[i].row+avg_scale/2>line_positions[5+j*10] && symbols[i].row+avg_scale/2<line_positions[5+j*10]+avg_scale*3/4)||(symbols[i].row+avg_scale/2>line_positions[9+j*10]-avg_scale/3 && symbols[i].row+avg_scale/2<line_positions[9+j*10]+avg_scale/3) ){
									symbols[i].pitch='G';
									break;
								}
								else if((symbols[i].row+avg_scale/2>line_positions[0+j*10]-avg_scale/3 && symbols[i].row+avg_scale/2 <=line_positions[0+j*10]+avg_scale/3)||(symbols[i].row+avg_scale/2>line_positions[3+j*10]+avg_scale/2 && symbols[i].row+avg_scale/2 <=line_positions[3+j*10]+avg_scale*3/4) ||(symbols[i].row+avg_scale/2>line_positions[6+j*10]-avg_scale/3 && symbols[i].row+avg_scale/2 <=line_positions[6+j*10]+avg_scale/3)||(symbols[i].row+avg_scale/2>line_positions[9+j*10]+avg_scale/2 && symbols[i].row+avg_scale/2 <=line_positions[9+j*10]+avg_scale*3/4)){
									symbols[i].pitch='F';
									break;
								}
								else if((symbols[i].row+avg_scale/2>line_positions[0+j*10]+avg_scale/2 && symbols[i].row+avg_scale/2 <=line_positions[0+j*10]+avg_scale*3/4)||(symbols[i].row+avg_scale/2>line_positions[4+j*10]-avg_scale/3 && symbols[i].row+avg_scale/2 <=line_positions[4+j*10]+avg_scale/3)||(symbols[i].row+avg_scale/2>line_positions[6+j*10]+avg_scale/2 && symbols[i].row+avg_scale/2 <=line_positions[6+j*10]+avg_scale*3/4)||((symbols[i].row+avg_scale/2>line_positions[9+j*10]+avg_scale*3/4 && symbols[i].row+avg_scale/2 <=line_positions[9+j*10]+avg_scale+1))){
									symbols[i].pitch='E';
									break;
								}
								else if((symbols[i].row+avg_scale/2>line_positions[1+j*10]-avg_scale/3 && symbols[i].row+avg_scale/2 <=line_positions[1+j*10]+avg_scale/3)||(symbols[i].row+avg_scale/2>=line_positions[4+j*10] && symbols[i].row+avg_scale/2 <=line_positions[4+j*10]+avg_scale*3/4)||(symbols[i].row+avg_scale/2>line_positions[7+j*10]-avg_scale/3 && symbols[i].row+avg_scale/2 <=line_positions[7+j*10]+avg_scale/3)){
										symbols[i].pitch='D';
										break;
									}
								else if((symbols[i].row+avg_scale/2>line_positions[1+j*10]+avg_scale/2 && symbols[i].row+avg_scale/2 <=line_positions[1+j*10]+avg_scale*3/4)||(symbols[i].row+avg_scale/2>line_positions[4+j*10]+avg_scale*3/4 && symbols[i].row+avg_scale/2 <=line_positions[4+j*10]+avg_scale)||(symbols[i].row+avg_scale/2>line_positions[5+j*10]-14 && symbols[i].row+avg_scale/2 <=line_positions[5+j*10]-avg_scale*3/4)||(symbols[i].row+avg_scale/2>line_positions[7+j*10]+avg_scale/2 && symbols[i].row+avg_scale/2 <=line_positions[7+j*10]+avg_scale*3/4)){
											symbols[i].pitch='C';
											break;
										}
								else if((symbols[i].row+avg_scale/2>line_positions[2+j*10]-avg_scale/3 && symbols[i].row+avg_scale/2 <=line_positions[2+j*10]+avg_scale/3)||(symbols[i].row+avg_scale/2>line_positions[4+j*10]+avg_scale && symbols[i].row+avg_scale/2 <line_positions[5+j*10]-avg_scale)||(symbols[i].row+avg_scale/2>line_positions[5+j*10]-9 && symbols[i].row+avg_scale/2 <line_positions[5+j*10]-avg_scale/3)||(symbols[i].row+avg_scale/2>line_positions[8+j*10]-avg_scale/3 && symbols[i].row+avg_scale/2 <line_positions[8+j*10]+avg_scale/3)){
												symbols[i].pitch='B';
												break;
											}
								else if((symbols[i].row+avg_scale/2>line_positions[2+j*10] && symbols[i].row+avg_scale/2 <=line_positions[2+j*10]+avg_scale*3/4) || (symbols[i].row+avg_scale/2>line_positions[5+j*10]-avg_scale/3 && symbols[i].row+avg_scale/2 <=line_positions[5+j*10]+avg_scale/3)||(symbols[i].row+avg_scale/2>line_positions[8+j*10]+avg_scale/2 && symbols[i].row+avg_scale/2 <=line_positions[8+j*10]+avg_scale*3/4)){
													symbols[i].pitch='A';
													break;
												}

			}

	}

}
}

SDoublePlane gaussianBlurkernel(float sigma,int filtersize)
{
  SDoublePlane kernel(filtersize,filtersize);
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


SDoublePlane sobel_edge_detector(const SDoublePlane &input,const SDoublePlane &xfilter,const SDoublePlane &yfilter )
{
  SDoublePlane output(input.rows(), input.cols());
  SDoublePlane x_output(input.rows(), input.cols());
  SDoublePlane y_output(input.rows(), input.cols());
  x_output=convolve_general(input,xfilter);
  y_output=convolve_general(input,yfilter);
  for(int i=0;i<output.rows();i++) {
    for(int j=0;j<output.cols();j++) {
    	if (!i || !j || i == output.rows() - 1 || j == output.cols() - 1) {
    		output[i][j] = 255;
    	} else {
    		output[i][j]=sqrt(pow (x_output[i][j],2)+pow (y_output[i][j],2));
    	}
    }
  }
  return output;
}

// Apply an edge detector to an image, returns the binary edge map
//
SDoublePlane find_edges(const SDoublePlane &input, double thresh){
  SDoublePlane output(input.rows(), input.cols());
  for(int i=0;i<input.rows();i++){
    for(int j=0;j<input.cols();j++){
      if(input[i][j]>thresh)
      	input[i][j]=255;
      else
        input[i][j]=0;
      output[i][j]=255-(unsigned char)(input[i][j]);
    }
  }


  // Implement an edge detector of your choice, e.g.
  // use your sobel gradient operator to compute the gradient magnitude and threshold

  return output;
}


class HoughLine{
  public:
  int votes, row, scale;
};

vector<HoughLine> houghtransform(const SDoublePlane &input_unfiltered){

  int filtersize=2;
  SDoublePlane row_filter(1,filtersize);
  for(int i=0; i<filtersize; i++)
  row_filter[0][i] = 1/((double) filtersize);
  SDoublePlane col_filter(filtersize,1);
  for(int j=0; j<filtersize; j++)
  col_filter[j][0] = 1/((double) filtersize);
  SDoublePlane input = convolve_separable(input_unfiltered, row_filter,col_filter);

  vector<HoughLine> res;
  SDoublePlane voting_space(input.rows(), input.cols()/4);
  int black = 150;
  int white = 200;
  int err = 5;
  int voteerr = 10;

  for(int r=0; r<input.rows(); r++){
    for(int c=0; c<input.cols(); c++){
      int lcol = c-1;
      int rcol = c+1;
      if( lcol<0 || rcol==input.cols()){
        continue;
      }
      if (input[r][c] <= black && input[r][rcol] <= black && input[r][lcol] <= black){ //check the pixel to the right and left to make sure we're on a line
        int four_below[4];
        int i = 0;
        int new_r = r+1;
        bool target_black = true;
        while(new_r < input.rows() && i < 4){
          if(target_black && input[new_r][c] <= black){
            four_below[i] = new_r;
            i++;
            target_black = false;
          }
          else if(!target_black && input[new_r][c] >= white){
            target_black = true;
          }
          new_r++;
        }
        if(i != 4){ //we didn't find 4 black pixels below the current row candidate
          continue; //move on to the next pixel
        }
        int diff = four_below[0] - r;
        int new_diff = four_below[1] - four_below[0];
        i=2;
        while(abs(new_diff - diff) <= err && i < 4){
          diff=new_diff;
          new_diff = four_below[i] - four_below[i-1];
          i++;
        }
        if(i==4){ //this pixel will vote for this row, with a spacing value of average diff
          int avg_diff = (four_below[3] - r)/4;
          voting_space[r][avg_diff]++;
        }

      }
    }
  }
  int sum = 0;
  int numrows = 0;
  for(int r=0; r<voting_space.rows(); r++){

    double max = *max_element(voting_space[r], voting_space[r]+voting_space.cols());
    double pos = max_element(voting_space[r], voting_space[r]+voting_space.cols()) - voting_space[r];
    if(max!=0){
      numrows++;
      sum+=max;
      HoughLine h;
      h.row = r;
      h.votes = (int)max;
      h.scale = (int)pos;
      res.push_back(h);
    }
  }
  float avg = ((float)sum)/numrows;

  vector<HoughLine> new_res;
  for(int i=0; i<res.size(); i++){
    if(res[i].votes>=avg*2){
      new_res.push_back(res[i]);
    }
  }

  res=new_res;
  new_res.clear();

  //find the peaks in voting space
  if(res.size()>6){
    for(int i=0; i<res.size(); i++){
      int votes = res[i].votes;
      int left = i-1<0 ? 0 : res[i-1].votes;
      int right = i+1==res.size() ? 0 : res[i+1].votes;
      if(abs(votes-left)>=voteerr && abs(votes-right)>=voteerr){ //a peak is greater than it's neighbor on either side
        new_res.push_back(res[i]);
      }
    }
  }

  res=new_res;
  new_res.clear();

  //ignore staves until after the current one has concluded (duplicate finds within a staff)
  int ignore_until=res[0].votes>res[1].votes?res[0].row:res[1].row;
  for(int i=0; i<res.size(); i++){
    if(res[i].row >= ignore_until){
      new_res.push_back(res[i]);
      ignore_until= res[i].row + res[i].scale * 5;
    }
  }
  res=new_res;

  // for(int i = 0; i< res.size(); i++){
  // 	cout<<"r "<<res[i].row<<": "<<res[i].votes <<" sp:"<<res[i].scale<<endl;
  // }
  return res;

}



SDoublePlane output_edge_detector(const SDoublePlane &input)
{
  SDoublePlane output(input.rows(), input.cols());
  SDoublePlane gaussian_kernel(21,21);
  SDoublePlane sobel_x_filter(3,3);
  SDoublePlane sobel_y_filter(3,3);
  SDoublePlane gaussianBlurImage(input.rows(),input.cols());
  float sigma=0.05;// Sigma can be - 0.01, 1.0, 5.0, and 10.0
  double thresh = 0;//Could be set to 200, if not dynamic
  int x_filter[3][3]={
    {-1,0,1,},
    {-2,0,2,},
    {-1,0,1}
  };
  int y_filter[3][3]={
    {1,2,1,},
    {0,0,0,},
    {-1,-2,-1}
  };
  for(int i=0; i<sobel_x_filter.rows(); i++){
    for(int j=0; j<sobel_x_filter.cols(); j++){
      sobel_x_filter[i][j] = x_filter[i][j];
    }
  }
  for(int i=0; i<sobel_y_filter.rows(); i++){
    for(int j=0; j<sobel_y_filter.cols(); j++){
      sobel_y_filter[i][j] = y_filter[i][j];
    }
  }
  gaussian_kernel=gaussianBlurkernel(sigma,21);
  gaussianBlurImage=convolve_general(input,gaussian_kernel);
  output=sobel_edge_detector(input,sobel_x_filter,sobel_y_filter);
  for (int i = 0; i < output.rows(); i++){
    for (int j = 0; j < output.cols(); j++){
      if (output[i][j] >= thresh){
        thresh = output[i][j];
      }
    }
  }

  output=find_edges(output,thresh/2.5);
  return output;

}

SDoublePlane binary_image(const SDoublePlane &input)
{
  SDoublePlane output(input.rows(), input.cols());

  for(int i=0;i<input.rows();i++){
    for(int j=0;j<input.cols();j++){
      if(input[i][j] <= 150) {
        output[i][j]=255;
      }
      else {
        output[i][j]=0;
      }
    }
  }

  return output;
}

double euclid(int x1, int y1, int x2, int y2) {
	return sqrt(pow(x1 - x2 ,2) + pow(y1 - y2 ,2));
}

void fill_edge_dp(vector<int> *nls, int r, int c, int originr, int originc, SDoublePlane &dp) {
	if (r == -1 || r == dp.rows() || c == -1 || c == dp.cols()) {
		return;
	}
	double dist = euclid(r,c,originr,originc);
	if (dist >= dp[r][c]) {
		return;
	}

	dp[r][c] = dist;
	nls->push_back(r);
	nls->push_back(c);
	nls->push_back(originr);
	nls->push_back(originc);
}

void expand_dp_edges(SDoublePlane &dp) {
	vector<int> ls;
	for(int i=0;i<dp.rows();i++) {
    for(int j=0;j<dp.cols();j++) {
      if(!dp[i][j]) {
    		ls.push_back(i);
    		ls.push_back(j);
    		ls.push_back(i);
    		ls.push_back(j);
      }
    }
  }
  while(ls.size() != 0) {
  	vector<int> nls;
  	for (int lsi = 0; lsi < ls.size(); lsi += 4) {
  		int r = ls[lsi];
  		int c = ls[lsi + 1];
  		int originr = ls[lsi + 2];
  		int originc = ls[lsi + 3];
  		for (int i = -1; i <= 1; i++) {
				for (int j = -1; j <= 1; j++) {
					if (i == 0 && j == 0) {
						continue;
					}
					fill_edge_dp(&nls, r + i, c + j, originr, originc, dp);
				}
			}
  	}
  	ls = nls;
  }

}

SDoublePlane closest_edge_pixel(const SDoublePlane &input_binary_image) {
  SDoublePlane dp(input_binary_image.rows(), input_binary_image.cols());
  for(int i=0;i<input_binary_image.rows();i++) {
    for(int j=0;j<input_binary_image.cols();j++) {
    	// DoublePlane is initialized to 0, set to very high value if not an edge.
      if (input_binary_image[i][j]) {
      	dp[i][j] = 0;
      } else {
      	dp[i][j]= 999999999;
      }
    }
  }
  expand_dp_edges(dp);

  return dp;
}

SDoublePlane FIT(const SDoublePlane &input_image,const SDoublePlane &template_image,vector<DetectedSymbol> &symbols,Type type)
{

  SDoublePlane output(input_image.rows(), input_image.cols());
  // SDoublePlane closest_edge_output(input_image.rows(), input_image.cols());
  SDoublePlane input_binary_image=binary_image(input_image);
  SDoublePlane closest_edge_output = closest_edge_pixel(input_binary_image);

  SDoublePlane template_binary_image=binary_image(template_image);

  double sums = 0;
  double count = 0;
  double min_match = 999999999;
  for(int i=0;i<closest_edge_output.rows();i++) {
    for(int j=0;j<closest_edge_output.cols();j++) {
    	if (i >= closest_edge_output.rows() - template_binary_image.rows() || j >= closest_edge_output.cols() - template_binary_image.cols()) {
    		output[i][j] = 999999999;
    		continue;
    	}
      double template_sum = 0;
      for(int k = 0; k < template_binary_image.rows(); k++) {
        for(int l = 0; l < template_binary_image.cols(); l++) {
        	if (template_binary_image[k][l]) {
        		template_sum +=   closest_edge_output[i+k][j+l];
        	}

        }
      }
      output[i][j] = template_sum;
      sums += template_sum;
      count++;
      if (template_sum < min_match) {
      	min_match = template_sum;
      }
    }
  }

  double avg = sums / count;

	double match_threshold = (min_match + 5)* 1.25;

  for(int i=0;i<output.rows();i++){
    for(int j=0;j<output.cols();j++){
      if(output[i][j] < match_threshold){
        DetectedSymbol s;
        s.row = i;
        s.col = j;
        s.width = template_binary_image.cols();
        s.height = template_binary_image.rows();
        s.type = type;
        symbols.push_back(s);
      }
    }
  }

  return output;
}

cv::Mat doubleplane_to_mat(const SDoublePlane &source){
  cv::Mat res(source.rows(), source.cols(), cv::DataType<float>::type);
  for(int r=0; r<source.rows(); r++){
    for(int c=0; c<source.cols(); c++){
      res.at<float>(r,c) = source[r][c];
    }
  }
  return res;
}

SDoublePlane mat_to_doubleplane(cv::Mat source){
  SDoublePlane res(source.rows, source.cols);
  for(int r=0; r<source.rows; r++){
    for(int c=0; c<source.cols; c++){
      res[r][c] = source.at<float>(r,c);
    }
  }
  return res;
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

  //Define seperable filter
  int filtersize=2;
  SDoublePlane row_filter(1,filtersize);
  for(int i=0; i<filtersize; i++)
  row_filter[0][i] = 1/((double) filtersize);
  SDoublePlane col_filter(filtersize,1);
  for(int j=0; j<filtersize; j++)
  col_filter[j][0] = 1/((double) filtersize);

  SDoublePlane gaussian = gaussianBlurkernel(1, 21);


  /**************
  ** File I/O **
  ***************/
  string input_filename(argv[1]);
  SDoublePlane input_image= SImageIO::read_png_file(input_filename.c_str());
  string template1_filename = "template1.png";
  SDoublePlane template1_image= SImageIO::read_png_file(template1_filename.c_str());
  string template2_filename = "template2.png";
  SDoublePlane template2_image= SImageIO::read_png_file(template2_filename.c_str());
  string template3_filename= "template3.png";
  SDoublePlane template3_image= SImageIO::read_png_file(template3_filename.c_str());

  SDoublePlane conv_template1_image = convolve_general(template1_image, gaussian);
  SDoublePlane conv_template2_image = convolve_general(template2_image, gaussian);
  SDoublePlane conv_template3_image = convolve_general(template3_image, gaussian);

  printf("File I/O completed\n");

  /**************
  * Question 2 *
  ***************/

  SDoublePlane output_1 = convolve_general(input_image, gaussian);
  printf("Mean filtering completed\n");

  /**************
  * Question 3 *
  ***************/

  SDoublePlane output_2 = convolve_separable(input_image,row_filter,col_filter);
  printf("Seperable convolution completed\n");

  /**************
  * Question 4 *
  ***************/

  vector<DetectedSymbol> symbols;
  write_detection_image("scores4.png", symbols, output_1);

  detect_symbols(output_1, conv_template1_image, symbols, NOTEHEAD);
  detect_symbols(output_1, conv_template2_image, symbols, QUARTERREST);
  detect_symbols(output_1, conv_template3_image, symbols, EIGHTHREST);
  vector<int> line_positions;
  if(input_filename!="music3.png"){ // NOTE: For music3, the detect_lines code written for Q4 doesn't return more than 10 lines. So in compute_pitch logic it results in segmentation fault.
		line_positions=detect_lines(input_image);
		line_positions.erase(line_positions.begin());
		line_positions.erase(line_positions.begin()+10);
		compute_pitch(line_positions,symbols,11);
  }
  write_detection_image("detected4.png", symbols, input_image);

  printf("Question 4 completed\n");

  /**************
  * Question 5 *
  ***************/

  SDoublePlane edge_input_image(input_image.rows(),input_image.cols());
  SDoublePlane edge_template_image1(template1_image.rows(),template1_image.cols());
  SDoublePlane edge_template_image2(template2_image.rows(),template2_image.cols());
  SDoublePlane edge_template_image3(template3_image.rows(),template3_image.cols());
  edge_input_image=output_edge_detector(input_image);
  edge_template_image1=output_edge_detector(template1_image);
  edge_template_image2=output_edge_detector(template2_image);
  edge_template_image3=output_edge_detector(template3_image);

  symbols.clear();
  write_detection_image("edges.png", symbols, edge_input_image);

  FIT(edge_input_image,edge_template_image1,symbols,NOTEHEAD);
  FIT(edge_input_image,edge_template_image2,symbols,QUARTERREST);
  FIT(edge_input_image,edge_template_image3,symbols,EIGHTHREST);

  write_detection_image("detected5.png", symbols, input_image);

  printf("Question 5 completed\n");

  /**************
  * Question 6 *
  ***************/

  vector<HoughLine> detected_lines = houghtransform(input_image);

  line_positions.clear();
  symbols.clear();
  //create symbol versions of the detected lines
  for(int i=0; i<detected_lines.size(); i++){
    int row = detected_lines[i].row;
    for(int j=0; j<5; j++){
      DetectedSymbol s;
      s.row =	row;
      s.col = 0;
      s.width = input_image.cols();
      s.height = 1;
      s.type = EIGHTHREST;
      symbols.push_back(s);
      line_positions.push_back(row);
      row+=detected_lines[i].scale;
    }
  }
  write_detection_image("staves.png", symbols, input_image);

  printf("Question 6 completed\n");

  /**************
  * Question 7 *
  ***************/

  int sum=0;
  for(int i=0; i<detected_lines.size(); i++){
    sum+=detected_lines[i].scale;
  }
  int avg_scale = sum/detected_lines.size();

  //Scale the template images before detection
  cv::Mat source;
  cv::Mat dest;

  cv::Size dest_size((int)(((float)avg_scale)/template1_image.rows() * template1_image.cols()), avg_scale);  //new size for notehead
  dest = cv::Mat(dest_size, cv::DataType<float>::type);
  source = doubleplane_to_mat(template1_image);
  cv::resize(source,dest,dest_size);
  SDoublePlane resize_template1_image = mat_to_doubleplane(dest);


  dest_size = cv::Size((int)(((float)avg_scale)/template1_image.rows() * template2_image.cols()), avg_scale*3);  //new size for quarterrest
  dest = cv::Mat(dest_size, cv::DataType<float>::type);
  source = doubleplane_to_mat(template2_image);
  cv::resize(source,dest,dest_size);
  SDoublePlane resize_template2_image = mat_to_doubleplane(dest);

  dest_size = cv::Size((int)(((float)avg_scale)/template1_image.rows() * template3_image.cols()), avg_scale*2);  //new size for eighthrest
  dest = cv::Mat(dest_size, cv::DataType<float>::type);
  source = doubleplane_to_mat(template3_image);
  cv::resize(source,dest,dest_size);
  SDoublePlane resize_template3_image = mat_to_doubleplane(dest);

  printf("Scaled template images\n");

  symbols.clear();

  detect_symbols(output_1,convolve_general(resize_template1_image,gaussian),symbols,NOTEHEAD);
  detect_symbols(output_1,convolve_general(resize_template2_image,gaussian),symbols,QUARTERREST);
  detect_symbols(output_1,convolve_general(resize_template3_image,gaussian),symbols,EIGHTHREST);
  compute_pitch(line_positions,symbols,avg_scale);
  write_detection_image("detected7.png", symbols, input_image);
  write_detection_txt("detected7.txt", symbols);

  printf("Question 7 completed\n");

}
