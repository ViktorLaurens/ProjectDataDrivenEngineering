# ProjectDataDrivenEngineering
In this repository, code can be found in relation to a project in which data driven analysis techniques are used on datasets containing data about the motion of various pendula. The objective of this project is to explore the vast methods and strategies that data and optimization offer. For this project, a (simple) mechanical system is chosen as the subject of the analysis, namely: a single pendulum with viscous friction. The pendulum serves as a classic example in physics that allows us to explore and understand fundamental concepts such as gravity, oscillatory motion, and conservation of energy. Despite its simplicity, the pendulum exhibits a rich variety of behaviors. This subject was chosen because it is intuitive, possible to build ourselves and, under certain assumptions, analytic solutions to the problem exist. This allows us to link the results from the data-driven methods to the theory and to explore the advantages and limitations of both and how these relate to each other.

## Table of Contents
- [Installation and usage](#installation)
- [Segments](#segments)
- [Documentation](#documentation)
- [License](#license)

## Installation and usage

All the code has been written in Google Colab so it can be simply copy pasted and run in a .ipynb file in Colab. 

## Segments

This repository provides the following segments:

### 1. Data-acquisition and Outlier Detection
The first step is to acquire the data. An experimental setup is built and the data is collected by tracking the ball using computer vision. However, performing these experiments are very time-consuming. A virtual data set is generated using the theoretical formulation of the problem. By using simulation to obtain
a diverse data set, we can efficiently explore a wide range of parameter variations and observe the corresponding effects on the pendulum’s motion. This approach saves significant time and resources compared to conducting physical experiments for each parameter combination.

### 2. PCA for dimensionality reduction
Principal Component Analysis (PCA) is a widely employed unsupervised learning technique in data analysis, particularly for large multidimensional data sets with numerous features. Its primary purpose is dimensionality reduction and feature extraction by transforming the correlated features into a set of uncorrelated features. Effective dimensionality reduction is achieved by transforming the data in the lower-dimensional reduced PCA space spanned by the first few principal components that sufficiently account for the variance in the data.

### 3. Linear regression for pendulum motion prediction
Linear regression provides a clear and intuitive understanding of how changes in independent variables influence a dependent variable. The coefficients of the linear regression model represent the direction and magnitude of these effects. This interpretability makes linear regression a valuable tool for gaining insights and making informed decisions. In addition, linear regression allows us to make predictions based on observed data. Once the model is trained, it can be used to estimate the values of the dependent variable for new data points. This predictive capability is valuable in various applications and will be used here to estimate the future position of the pendulum based on previous positions. 

### 4. Modal analysis and proper orthogonal decomposition
Modal analysis is another possible way to analyze data. The main idea of Modal analysis is to decompose a data set as a linear combination of its modes. Those modes are a linear combination of elementary contributions. Principally, it is also linked with PCA. The main difference between PCA and modal analysis or rather proper orthogonal decomposition (POD), lies in the approach. POD a via singular value decomposition (SVD) is considered in the literature as a method of POD as well as PCA and are also seen as equivalent to each other. The relevant difference for this report lies in the implementation. While PCA focuses on the principal components, POD via SVD focuses on the spatial and temporal distribution of the data.

### 5. Classification and regression using a neural network
Today, neural networks (NNs) are a hot topic. Their popularity has increased tremendously over recent years thanks to an increase in computational power and the availability of large data sets. They are a subset of supervised learning where the algorithm is trained on labeled input data. Once trained, the network can be used for recognizing patterns, classifying data, and making predictions. A NN with sufficient hidden neurons and appropriate activation functions is considered to be a universal function approximation. It can approximate any arbitrary function. However, tuning the NN, for example, choosing the number of hidden layers and neurons per layer, is still challenging. Generally, the more complex the problem, the more layers are required.

## Documentation

Comprehensive documentation about the different methods used and the findings and insights obtained by using them can be found here:
[MA1_Project__Introducing_Adaptive_Extension_to_Bidirectional_Rapidly_Exploring__Random_Trees_for_Multi_Robot_Systems.pdf](https://github.com/ViktorLaurens/MA1_Project/files/11470170/MA1_Project__Introducing_Adaptive_Extension_to_Bidirectional_Rapidly_Exploring__Random_Trees_for_Multi_Robot_Systems.pdf)


## License

This project is licensed under the [MIT License](LICENSE).

The MIT License is a permissive open-source license that allows you to use, modify, and distribute this software for any purpose, both commercially and non-commercially. The license text can be found in the [LICENSE](LICENSE) file.

### Permissions
- Commercial use: ✔️
- Modification: ✔️
- Distribution: ✔️
- Private use: ✔️

### Limitations
- Liability: ❌
- Warranty: ❌

By using this software, you agree to the terms and conditions of the MIT License.

For more information about the MIT License, please visit [opensource.org/licenses/MIT](https://opensource.org/licenses/MIT).
