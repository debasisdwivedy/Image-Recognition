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

    int c=col_filter.rows()/2;

    int imageRows=input.rows();
    int imageCols=input.cols();

    for(int i=c;i<imageRows-c;i++){
  	  for(int j=c;j<imageCols-c;j++){
  		  double temp=0;

  		  	  for(int v=-c;v<c;v++){
  		  			temp=temp+row_filter[0][v+c] * input[i][j-v];

  		  		  }

  			output_temp[i][j]=temp;
  	  }

    }

    for(int i=c;i<imageRows-c;i++){
      	  for(int j=c;j<imageCols-c;j++){
      		  double temp=0;

      		  	  for(int u=-c;u<c;u++){
      		      		 temp=temp+col_filter[u+c][0] * output_temp[i-u][j];

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

  for(int i=c;i<imageRows-c;i++){
	  for(int j=c;j<imageCols-c;j++){
		  double temp=0;
		  for(int u=-c;u<c;u++){

  			  for(int v=-c;v<c;v++){
  				  temp=temp+filter[u+c][v+c] * input[i+u][j+v];

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
			if(image[i][j]>70){
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
		cout<< "row: " <<i <<endl;
		for(int j=0;j<image.cols();j++){
			cout<< image[i][j]<<" ";
		}
		cout<< endl <<endl;
	}
	cout<<"display finished"<<endl;
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
	for(int i=d;i<binarized_input_image.rows()-d;i++){
		for(int j=c;j<binarized_input_image.cols()-c;j++){
			temp=0;
			for(int k=0;k<binarized_template_image.rows();k++){
				for(int l=0;l<binarized_template_image.cols();l++){
					temp=temp+binarized_input_image[i+k-d][j+l-c] * binarized_template_image[k][l]+(1-binarized_input_image[i+k-d][j+l-c]) * (1-binarized_template_image[k][l]);
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
				      s.row =i-binarized_template_image.rows()/2-4;
				      s.col =j-binarized_template_image.cols()/2-2;
				      s.width = binarized_template_image.cols()+4;
				      s.height = binarized_template_image.rows()+6;
				      s.type = type;
				      symbols.push_back(s);
				      block_image_area(hamming_matrix,s.row,s.col,s.width,s.height);
			}
		}
	}




}
/*
 * This method detect the co-ordinates of the staff lines in the image.
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
		for (int i = 1; i < input.rows() - 1; i++) {
			for (int j = 1; j < input.cols() - 1; j++) {
				int temp = 0;
				for (int u = 0; u < 3; u++) {
					for (int v = 0; v < 3; v++) {
						temp = temp
								+ input[i + u - 1][j + v - 1]
										* horizontal_line_kernel[u][v];
					}
				}
				if (temp<-460)
					output_lines[i][j] = 0;
				else
					output_lines[i][j] = 255;
			}
		}


//	return output_lines;
	vector<int> indices;
	  for(int i=0;i<output_lines.rows();i++){
		  int count_black=0;
		  for(int j=0;j<output_lines.cols();j++){
			  if(output_lines[i][j]>=0 && output_lines[i][j]<=33){
				  count_black++;
			  }
		  }
		  if(count_black>0.20*output_lines.cols()){
//
			  indices.push_back(i);

			  }
//		  else{
//			  for(int j=0;j<output_lines.cols();j++){
//			  			  output_lines[i][j]=255;
//			  }
//		  }

	  }
//	 return output_lines;
return indices;

}

/*
 * This method computes the pitch of the detected notes.
 */
void compute_pitch(vector<int> line_positions,vector<DetectedSymbol> &symbols){
for(int i=0;i<symbols.size();i++){
	for(int j=0;j<line_positions.size()/10-1;j++){
		if(symbols[i].row+9 <line_positions[1+j*10] || (symbols[i].row+9>line_positions[4+j*10]-3 && symbols[i].row+9<line_positions[4+j*10]+3) ||(symbols[i].row+9>line_positions[6+j*10]+5 && symbols[i].row+9<line_positions[6+j*10]+8)||(symbols[i].row+9>line_positions[10+j*10]-3 && symbols[i].row+9<line_positions[10+j*10]+3) ){
				symbols[i].pitch='G';
			}
			else if((symbols[i].row+9>line_positions[1+j*10]-3 && symbols[i].row+9 <=line_positions[1+j*10]+3)||(symbols[i].row+9>line_positions[4+j*10]+5 && symbols[i].row+9 <=line_positions[4+j*10]+8) ||(symbols[i].row+9>line_positions[7+j*10]-3 && symbols[i].row+9 <=line_positions[7+j*10]+3)||(symbols[i].row+9>line_positions[10+j*10]+5 && symbols[i].row+9 <=line_positions[10+j*10]+8)){
				symbols[i].pitch='F';
			}
			else if((symbols[i].row+9>line_positions[1+j*10]+5&& symbols[i].row+9 <=line_positions[1+j*10]+8)||(symbols[i].row+9>line_positions[5+j*10]-3 && symbols[i].row+9 <=line_positions[5+j*10]+3)||(symbols[i].row+9>line_positions[7+j*10]+5 && symbols[i].row+9 <=line_positions[7+j*10]+8)||((symbols[i].row+9>line_positions[10+j*10]+8 && symbols[i].row+9 <=line_positions[10+j*10]+14))){
				symbols[i].pitch='E';
			}
			else if((symbols[i].row+9>line_positions[2+j*10]-3 && symbols[i].row+9 <=line_positions[2+j*10]+3)||(symbols[i].row+9>line_positions[5+j*10]+5&& symbols[i].row+9 <=line_positions[5+j*10]+8)||(symbols[i].row+9>line_positions[8+j*10]-3&& symbols[i].row+9 <=line_positions[8+j*10]+3)){
					symbols[i].pitch='D';
				}
			else if((symbols[i].row+9>line_positions[2+j*10]+5 && symbols[i].row+9 <=line_positions[2+j*10]+8)||(symbols[i].row+9>line_positions[5+j*10]+8 && symbols[i].row+9 <=line_positions[5+j*10]+10)||(symbols[i].row+9>line_positions[6+j*10]-14 && symbols[i].row+9 <=line_positions[6+j*10]-8)||(symbols[i].row+9>line_positions[8+j*10]+5 && symbols[i].row+9 <=line_positions[8+j*10]+8)){
						symbols[i].pitch='C';
					}
			else if((symbols[i].row+9>line_positions[3+j*10]-3 && symbols[i].row+9 <=line_positions[3+j*10]+3)||(symbols[i].row+9>line_positions[5+j*10]+16 && symbols[i].row+9 <line_positions[5+j*10]+20)||(symbols[i].row+9>line_positions[6+j*10]-9 && symbols[i].row+9 <line_positions[6+j*10]-3)||(symbols[i].row+9>line_positions[9+j*10]-3 && symbols[i].row+9 <line_positions[9+j*10]+3)){
							symbols[i].pitch='B';
						}
			else if((symbols[i].row+9>line_positions[3+j*10]+5 && symbols[i].row+9 <=line_positions[3+j*10]+8) || (symbols[i].row+9>line_positions[6+j*10]-3 && symbols[i].row+9 <=line_positions[6+j*10]+3)||(symbols[i].row+9>line_positions[9+j*10]+5 && symbols[i].row+9 <=line_positions[9+j*10]+8)){
								symbols[i].pitch='A';
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

  printf("read input file\n");

  // test step 2 by applying mean filters to the input image
  int filtersize=2;
  SDoublePlane mean_filter(filtersize,filtersize);
  for(int i=0; i<filtersize; i++)
    for(int j=0; j<filtersize; j++)
      mean_filter[i][j] = 1/((double) filtersize*filtersize);
  SDoublePlane output_1 = convolve_general(input_image, mean_filter);

  printf("mean filters passed\n");
	
	  // Separable filter
  SDoublePlane row_filter(1,filtersize);
  for(int i=0; i<filtersize; i++)
	  row_filter[0][i] = 1/((double) filtersize);
  SDoublePlane col_filter(filtersize,1);
  for(int j=0; j<filtersize; j++)
  	  col_filter[j][0] = 1/((double) filtersize);
  SDoublePlane output_2 = convolve_separable(input_image, row_filter,col_filter);

  printf("seperable filter passed\n");

  //Hough lines
  vector<HoughLine> detected_lines = houghtransform(input_image);
  printf("Hough transform passed\n");

  

  string template1_filename = "template1.png";
  SDoublePlane template1_image= SImageIO::read_png_file(template1_filename.c_str());

  printf("opened template1.png\n");

  string template2_filename = "template2.png";
  SDoublePlane template2_image= SImageIO::read_png_file(template2_filename.c_str());

  printf("opened template2.png\n");

  string template3_filename = "template3.png";
  SDoublePlane template3_image= SImageIO::read_png_file(template3_filename.c_str());

  printf("opened template3.png\n");


  vector<DetectedSymbol> symbols;

  //add detected lines
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
	  	row+=detected_lines[i].scale;
  	}
  }

  detect_symbols(input_image,template1_image,symbols,NOTEHEAD);
  detect_symbols(input_image,template2_image,symbols,QUARTERREST);
	detect_symbols(input_image,template3_image,symbols,EIGHTHREST);

  vector<int> line_positions=detect_lines(input_image);
  compute_pitch(line_positions,symbols);






  write_detection_txt("detected.txt", symbols);
  write_detection_image("convolution_general.png", symbols, output_1);
  write_detection_image("convolution_seperable.png", symbols, output_2);
  write_detection_image("detected.png", symbols, input_image);


}
