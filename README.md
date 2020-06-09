# PCA_fruit_images

Principal Component Anaylsis on Fruit Image Data

- Given RGB images of 100 fruits, *Principal Component Anaylsis* was performed, and the closest representations 
was obtained as a linear combination of the mean and the four eigenvectors corresponding to the four most significant 
eigenvalues of the covariance matrix. A MultiVariate Gaussian was fitted on the entire dataset of about 100 images.
- New fruit images were generated by random sampling, using the closest representations, which were 
distinct from any fruit in the dataset, but representative of the dataset.

Instructions to run code and theory is in report, and all the graphs are in results

### Detailed Desription:

**Principal Component Analysis (PCA) for Fruit Images:**

Consider the dataset provided within the folder “data_fruit”

Each datum is an image of size 80 × 80 pixels with 3 color channels red (R), green (G), and blue
(B), i.e., a 80 × 80 × 3 array. For PCA, each image should be resized to a vector of length 19200 .

For visualization, reshape each vector back to a RGB image of size 80 × 80 pixels using the
function reshape(), followed by a shift and rescaling of the values into the range [0, 1] , followed by
displaying the matrix using the function image().

• Find the mean μ , the covariance matrix C , and the the first 4 eigenvectors of C. 
Display the mean and the eigenvectors as images (side by side, in the same figure);  
Find the first 10 eigenvalues, sort them, and plot their values on a graph. 

• For each fruit image in the dataset, finds its closest representation as a
linear combination of the mean and the first 4 eigenvectors. Use the measure of closeness
as the Frobenius norm of the difference. Describe the algorithm used to produce this closest
representation in mathematical terms and describe the logic behind your algorithm. Display the
original fruit image and its closest representation, as images (side by side, in the same figure).

• Using all of the top 4 eigenvectors and the mean image, sample random
images from the multivariate Gaussian to generate new images of "fruit". Display three such
images that are distinct from any image in the given the dataset, but are representative of the
dataset and can be considered as that of a new / generated fruit.

