# Load necessary libraries
library(gbm)
library(ROCR)
library(caTools)

# Function to calculate the AUC (Area Under Curve) for test data
GBMTestAUC = function(model, testdata, ActualResponse) {
  # Predict probabilities for test data using the GBM model
  VirusGBMTestPred = predict.gbm(model, newdata=testdata, n.trees=model$n.trees, type="response")
  
  # Generate prediction object for performance evaluation
  VirGBMPred = prediction(VirusGBMTestPred, ActualResponse)
  
  # Calculate performance metrics (True Positive Rate vs. False Positive Rate)
  VirGBMPerf = performance(VirGBMPred, measure = "tpr", x.measure = "fpr")
  
  # Calculate AUC (Area Under the Curve)
  VirGBMauc = performance(VirGBMPred, measure="auc")
  
  # Extract AUC value
  a = as.numeric(VirGBMauc@y.values)
  
  return(c(a)) # Return AUC score
}

# Function to build multiple GBM models and store results
ModelBuilder = function(data, treeDepth=3, Initial.Trees=4000, ModelIter=10, LearningRate=0.001) {
  Models = list()  # List to store models
  SelectedTrees = c()  # Stores selected number of trees per model
  
  trainAUC = c()  # Stores AUC for training data
  testAUC = c()  # Stores AUC for test data
  testAUCforCountriesWithVectors = c()  # Placeholder for additional AUC calculation (currently commented out)
  
  RelInfluence = list()  # Stores relative importance of features for each model
  PartialDependence = list()  # Stores partial dependence data for each feature
  PredictedResponse = matrix(nrow=nrow(data), ncol=ModelIter)  # Matrix to store predicted responses for each iteration
  
  # Main loop to train and evaluate models across iterations
  for(i in 1:ModelIter) {
    print(i)  # Print iteration number for tracking
    set.seed(i)  # Set seed for reproducibility
    
    # Split data into training (80%) and test (20%) sets
    split = sample.split(data[,"VirusPresent"], 0.8)
    TrainingSet = data[which(split),]
    TestSet = data[which(!split),]
    
    # Ensure column names match original data
    names(TrainingSet) = names(data)
    names(TestSet) = names(data)
    
    ### GBM Model Training
    # Train a GBM model with specified parameters
    GbmModel = gbm(VirusPresent ~ .,  # Predict VirusPresent
                   data = select(TrainingSet, -Virus, -Country),  # Exclude Virus and Country columns
                   interaction.depth = treeDepth,  # Depth of trees
                   n.trees = Initial.Trees,  # Initial number of trees
                   distribution = "bernoulli",  # Binary classification
                   cv.folds = 10,  # Cross-validation folds
                   shrinkage = LearningRate)  # Learning rate
    
    # Select optimal number of trees based on cross-validation
    perf = as.character(gbm.perf(GbmModel, plot.it = FALSE, method = "cv"))
    SelectedTrees[i] = GbmModel$n.trees  # Store selected tree count
    
    ## AUC Calculations
    # Calculate AUC for training and test data
    trainAUC[i] = GBMTestAUC(GbmModel, TrainingSet, TrainingSet$VirusPresent)
    testAUC[i] = GBMTestAUC(GbmModel, TestSet, TestSet$VirusPresent)
    
    ## Feature Importance Calculation
    # Calculate feature importance and store it in a list
    RelInf = summary.gbm(GbmModel, plotit = FALSE, order = F, method = permutation.test.gbm)
    RelInfluence[[i]] = RelInf[order(RelInf[,1]),]  # Order by feature names
    
    ## Partial Dependence Calculation
    # For each feature in the importance list, calculate partial dependence
    for(j in 1:nrow(RelInf)) {
      # Check if feature is numeric before calculating partial dependence
      if(is.numeric(TrainingSet[,as.character(RelInf[j,1])])) {
        # Get partial dependence values for the feature
        a = plot.gbm(GbmModel, i.var = as.character(RelInf[j,1]), return.grid = TRUE, type = "response")
        a$model = rep(i, nrow(a))  # Add model iteration info
        
        # Initialize or append partial dependence data for each feature
        if(i == 1) {
          PartialDependence[[j]] = a 
          names(PartialDependence[j]) = as.character(RelInf[j,1])
        } else {
          PartialDependence[[j]] = rbind(PartialDependence[[j]], a)  
        }
      }
    }
    
    # Assign names to PartialDependence list based on data column names
    names(PartialDependence) = data %>% select(-Virus, -Country, -VirusPresent) %>% names()
    
    ## Prediction for Complete Dataset
    # Generate predictions for all data and store in PredictedResponse matrix
    Pred = predict.gbm(GbmModel, data, type = "response")
    PredictedResponse[,i] = Pred
  }
  
  ### Aggregation of Relative Influence Across Models
  MeanInfluence = c()  # Mean relative importance
  SDinfluence = c()  # Standard deviation of importance
  RelativeInfluence = RelInfluence[[1]]  # Initialize with first model's importance
  
  # Loop to combine importance from all models
  for(k in 2:ModelIter) {
    RelativeInfluence = cbind(RelativeInfluence, RelInfluence[[i]][,2])
  }
  # Calculate mean and standard deviation for each feature
  for(l in 1:nrow(RelativeInfluence)) {
    SDinfluence[l] = sd(subset(RelativeInfluence, select = -1)[l,], na.rm = T)
  }
  MeanInfluence = rowMeans(subset(RelativeInfluence, select = -1), na.rm = T)
  
  # Create data frame with mean and SD of relative influence for each feature
  RelativeInfluence = data.frame("Variable" = RelativeInfluence[,1], MeanInfluence, SDinfluence)
  
  ### Aggregation of Predictions Across Models
  MeanResponse = c()  # Mean predictions for each row
  ResponseSD = c()  # Standard deviation of predictions
  Country = data$Country  # Store country info
  Virus = data$Virus  # Store virus info
  ActualResponse = data$VirusPresent  # Actual response
  
  # Calculate mean and SD of predictions for each row
  for(i in 1:nrow(PredictedResponse)) {
    MeanResponse[i] = mean(PredictedResponse[i,])
    ResponseSD[i] = sd(PredictedResponse[i,])
  }
  
  # Create final data frame with predictions, uncertainty, and actual responses
  PredictedResponseFrame = data.frame(Virus, Country, MeanResponse, ResponseSD, ActualResponse)[order(MeanResponse, decreasing = TRUE),]
  
  # Return list of results
  return(list("SelectedTrees" = SelectedTrees,
              "TrainingAUC" = trainAUC,
              "TestingAUC" = testAUC,
              "PartialDependence" = PartialDependence,
              "PredictedResponse" = PredictedResponseFrame,
              "RelativeInfluence" = RelativeInfluence))
}

# Example implementation of primary model using filtered dataset
PrimaryModel = ModelBuilder(data = filter(VirusCountryPairs,
                                          AnthroVectorsInCountry > 0 | NonAnthroVectorsInCountry > 0) %>%
                              select(-AnthroVectorsAndCarriersInCountry,                         
                                     -NonAnthroVectorsAndCarriersInCountry,
                                     -Log10.NCBI.Sequences),
                            ModelIter = 10)
