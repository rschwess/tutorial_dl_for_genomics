
# DeepLearning for Genomics Tutorial

In this tutorial we will use a convolutional neuronal network to address a fairly basic but common problem in genomics. Given a set of sequences belonging to different classes, what are the characteristics in the DNA sequence that let us distinguish the classes. For example, given a set of promoter and enhancer sequences, we could ask if there are any patterns in the DNA that let us distinguish between the two. 

For a 'simple' example, think of two ChIP-seq experiments for two different transcription factors. After analysing the ChIP-seq data, we know at what positions in the genome the two factors bind and the majority of binding sites might be distinct. If we extract the underlying sequences and search for DNA patterns that are enriched in the respective sets, we get an idea what DNA sequences the transcription factor might bind and/or which co-factors influence their binding.

Traditionally, people use motif discovery tools for finding overrepresented words and motifs. However, if we move to slightly more complicated questions these methods quickly reach their limits and machine learning approaches become more promising.

A more complicated question: If we have multiple sets of enhancers that are active in different tissues and cell types and we have the underlying sequences, can we figure out what sequence patterns are characteristic for what activity? And once we know that can we infer which factors are common and which are tissue specific?

------

For our test dataset we have a simplfied, simulated version of such a task. We simulated 40,000 DNA sequences of length 200 bp. We split them into 4 enhancer classes and populated them with transcription factor binding motifs and other DNA patters to make them distinguishable. However, some motifs are shared between classes, they may overlap each other and are not necessarily perfect matches to the text book motifs. Thats much more how regulatory DNA actually looks like :)!

We will use keras to build and train a small convolutional neuronal network to classify our enhancer sequences. Once this network is trained well, we can than investigate how the network has learned to distinguish between the classes and try to relate this back to transcriptions factor motifs and so on.

You can either run everything in an interactive python or ipython session or (especially later when optimizing) just adjust and run python dl_intro.py in the terminal.

-----

# Data Set Up and Import

Lets start by looking at the data. We have 40,000 sequences and they are all labeled with their respective class. We already split them up into taining, test and validation set. It is also a good idea to check if our classes are roughly equally distributed across our different sets.


```bash
%%bash

# check lines per set
echo "Numbers:"
wc -l ./data/pwm*

# check data format
echo -ne "\nFormat\n"
head -n 3 ./data/pwm_seq_200bp_test_set.txt

# check class representation
echo -ne "\nClass Representations:\n"

echo -ne "\nTraining:\n"
cut -f 1 ./data/pwm_seq_200bp_train_set.txt | sort | uniq -c
echo -ne "\nTest:\n"
cut -f 1 ./data/pwm_seq_200bp_test_set.txt | sort | uniq -c
echo -ne "\nValidation:\n"
cut -f 1 ./data/pwm_seq_200bp_valid_set.txt | sort | uniq -c
```

    Numbers:
       1000 ./data/pwm_seq_200bp_test_set.txt
      38000 ./data/pwm_seq_200bp_train_set.txt
       1000 ./data/pwm_seq_200bp_valid_set.txt
      40000 total
    
    Format
    3	ATGGCTGATAATGACGATTGTACAGATGGTGGATGAGATTGCCTCGTCCCGGCAGCATTACCCCCTGGTGGCAACGGCCACCAGGGGGCAATAAATCTGTGTCTTATCTCCGAGACCAAACAATTCCACAGCCTCTTATACAGCACCGAATGGACCGCCCCCTGGTGGCCAGGTATCGTCGAGGGCTCAATTAAACTCCT
    1	GCAGGCATTATGAGGTAATAAACTCAGCGCGTGTTGAGATAAGATTCTAAGCGGCGCGCGCGCGCGACCGCGAGAAGTGGAGATTAAGCGCGCTAATGGTGTGTCCGATAGTCACGTGTCCGCGCGGCGCGCGCCATGTATGTTCTGTTCTGCGCGCCGCGCTTTGCGCGCGCGCTTGGTATATAAAGCTGGGTTTTAAT
    1	GGCGCGCCTGGCATTTCTTAGAGAGGCGCGCAATACAACGAGAATCACCTAGAAGCCGTGTCTGTTGCTTATCACCGTTCGCCTAGGCCGCACGGGCACGTGGGTCTCCCGTTCCCTCAATCCTAACAGAAGCGCGCTAAGTCGTCGTTGGCTCTCTTACTAGCAGCGCGCCTGTACTAACCCGGCACTCGGCGGTGGGC
    
    Class Representations:
    
    Training:
       9489 0
       9513 1
       9508 2
       9490 3
    
    Test:
        268 0
        237 1
        243 2
        252 3
    
    Validation:
        243 0
        250 1
        249 2
        258 3


-----

Looks good, lets move straight in! We will use keras with tensorflow as its backend. Keras is ideal for quickly writing down and prototyping networks in just a few lines of code. The documentation site will be usefull throught the tutorial https://keras.io/.

We import keras and the relevant layers and operations we need. (Don't worry about the anoying warning if you run under python 3.6)


```python
import keras
from keras.models import Sequential
from keras.layers import Dense, Flatten
from keras.layers import Conv1D, MaxPooling1D
import numpy as np
import os
```

    Using TensorFlow backend.
    /home/ron/anaconda3/envs/dl_tutorial/lib/python3.6/importlib/_bootstrap.py:219: RuntimeWarning: compiletime version 3.5 of module 'tensorflow.python.framework.fast_tensor_util' does not match runtime version 3.6
      return f(*args, **kwds)


I wrote two helper functions to convert the sequences into hot encoded sequences and a wrapper to read in and assemble the data. Feel free to skip over this but you might want to have a quick look and understand how we format the data. The hot encoding transforms the sequence into an X x 4 array whith rows corresponding to the sequence position and the columns representing the 4 DNA bases. The respective base column that matches the sequence at that position is 1 the rest 0.


```python
# Helper Function  get hotcoded sequence
def get_hot_coded_seq(sequence):
    """Convert a 4 base letter sequence to 4-row x-cols hot coded sequence"""
    # initialise empty
    hotsequence = np.zeros((len(sequence),4))
    # set hot code 1 according to gathered sequence
    for i in range(len(sequence)):
        if sequence[i] == 'A':
            hotsequence[i,0] = 1
        elif sequence[i] == 'C':
            hotsequence[i,1] = 1
        elif sequence[i] == 'G':
            hotsequence[i,2] = 1
        elif sequence[i] == 'T':
            hotsequence[i,3] = 1
    # return the numpy array
    return hotsequence

# Helper function to read in the labels and seqs and store as hot encoded np array
def read_data(infile):
    # read file in
    with open(infile, "r") as f:
        seqs = []
        labels = []
        for i,l in enumerate(f):
            l = l.rstrip()
            l = l.split("\t")
            seqs.append(l[1])
            labels.append(l[0])
    # make labels np.array
    labels = np.array(labels)
    # convert to one_hot_labels
    hot_labels = keras.utils.to_categorical(labels, num_classes=4)
    # make seqs np.array
    hot_seqs = np.zeros( (len(seqs), 200, 4) )
    # fill with hot encoded sequences
    for j in range(len(seqs)):
        hotsequence = get_hot_coded_seq(seqs[j])
        hot_seqs[j,] = hotsequence
    return hot_labels, hot_seqs
```

Now we can read in the data.


```python
# read data --------------------------------------------------------------------
train_file = "./data/pwm_seq_200bp_train_set.txt"
train_labels, train_seqs = read_data(train_file)
valid_file = "./data/pwm_seq_200bp_valid_set.txt"
valid_labels, valid_seqs = read_data(valid_file)
```

Lets check how the data look like after we read and hot encoded it.


```python
# check shapes
print("Train Seq Shape", train_seqs.shape)
print("Train Label Shape", train_labels.shape)

# check data format
print("Labels Format:")
print(train_labels[1:5])
print("Seq Format (first 10 bp):")
print(train_seqs[1, 1:10,:])
```

    Train Seq Shape (38000, 200, 4)
    Train Label Shape (38000, 4)
    Labels Format:
    [[ 1.  0.  0.  0.]
     [ 0.  0.  1.  0.]
     [ 0.  1.  0.  0.]
     [ 0.  0.  0.  1.]]
    Seq Format (first 10 bp):
    [[ 0.  0.  0.  1.]
     [ 0.  0.  0.  1.]
     [ 0.  0.  1.  0.]
     [ 0.  1.  0.  0.]
     [ 0.  0.  1.  0.]
     [ 0.  1.  0.  0.]
     [ 0.  0.  1.  0.]
     [ 0.  1.  0.  0.]
     [ 0.  1.  0.  0.]]


The labels have a one hot encoding where every column represents a different class and in this case only one class can be active at a time.

The sequences have shape [sample x sequence_length x basis]. As a comparison a set of 2D images would have the dimensions [sample x pixel_rows x pixel_columns x colour_channels]. A grey scale picture would have only one channel while RGB images have three. We can thus think of our sequence as a 1D image with 4 channels.

----

# Building the Network

We now define our network. We first set some global and network architecture options and put them all together in the keras sequential mode. The sequential mode is an easy wrapper for linearly stacked networks that makes your code even more concise. We just define the model to be sequential and than add/stack layer after layer. Here we use a simple convolutional architecture. 

* Our first layer is a 1D convolution over the input:
   * We use a 1D convolution because we only want the filter to move along the sequence axis and map the channels to the hidden units.
   * We start with 10 hidden units or filters or kernels which are all of length 5 (bp)
   * We use the RELU activation function
   * We also define the input shape and how to pad the input if necessary (see doc.)
* We next perform max pooling where we take the maximum of a window of 5 consecutive activation values
    * This reduces the data dimension, thus simplifying the model and speeding up further computations
    * But it also enforces some extend of positional invariance into our model. For example, if we have a match to transcription factor motif in our sequence, we don't necessarily care where exactly this motif lies and a few bp up- or downstream shouldn't make a difference to our predictions.
* We then "Flatten" the activation values to a 1 dimensional vector
* And apply a fully connected or "Dense" layer connecting every value in the 1D vector to every class prediction
    * we use the sigmoid (softmax) activation function to perform effectively a multinomial logistic regression 

The number of hidden units, the size of the kernel, the pooling size, but also the number and types of layers we use in the network are usually called hyperparameters. The most work in DL usually comes down to finding the right hyperparameters that let our network training converge and that give us the best possible (at least the best we are able to find) accuracies.

We set some reasonable choices to begin with. But your task will be to play with these hyperparameters and see how well you can tune the model with a few adjustments.


```python
# global options
num_classes = 4

# network architecture options
conv1_hidden_units = 10
conv1_filter_size = 5
maxpool1_width = 5

# construct the model ----------------------------------------------------------
model = Sequential()
model.add(Conv1D(conv1_hidden_units, kernel_size=(conv1_filter_size), activation='relu', input_shape=(200, 4), padding='same'))
model.add(MaxPooling1D(pool_size=maxpool1_width))
model.add(Flatten())
model.add(Dense(num_classes, activation='softmax'))
```

Next we compile the model. We use adam as our optimizer. Since the classes are mutually exclusive we select the binary_crossentropy as our loss function and we want to monitor the accuracy during training.

We also print a summary of our network telling us the data shapes throught the network and sumarizing the number of trainable parameters in our model.


```python
# compile
model.compile(optimizer='adam',
              loss='binary_crossentropy',
              metrics=['accuracy'])

# print model summary
model.summary()
```

    _________________________________________________________________
    Layer (type)                 Output Shape              Param #   
    =================================================================
    conv1d_1 (Conv1D)            (None, 200, 10)           210       
    _________________________________________________________________
    max_pooling1d_1 (MaxPooling1 (None, 40, 10)            0         
    _________________________________________________________________
    flatten_1 (Flatten)          (None, 400)               0         
    _________________________________________________________________
    dense_1 (Dense)              (None, 4)                 1604      
    =================================================================
    Total params: 1,814
    Trainable params: 1,814
    Non-trainable params: 0
    _________________________________________________________________


Note that the fully connected layer has in comparison way more parameters to the convolutional layer.

-----

# Training

Now that we have our model or graph set up we can train it. We feed the model with our training sequences and labels, we define a batch size (since we are training in batch mode) and set the number of epochs (cycles through the training data) we want to train for. Five epochs should be fine for us feel free to ramp this up a bit and see if you get improvements or if the learning plateus quickly. 


```python
# Training Options
batch_size = 100
epochs = 5

# Train ------------------------------------------------------------------------
model.fit(train_seqs, train_labels,
          batch_size=batch_size,
          epochs=epochs,
          verbose=1,
          validation_data=(valid_seqs, valid_labels))
```

    Train on 38000 samples, validate on 1000 samples
    Epoch 1/5
    38000/38000 [==============================] - 4s 100us/step - loss: 0.3452 - acc: 0.8439 - val_loss: 0.2210 - val_acc: 0.9100
    Epoch 2/5
    38000/38000 [==============================] - 4s 100us/step - loss: 0.1913 - acc: 0.9220 - val_loss: 0.1671 - val_acc: 0.9315
    Epoch 3/5
    38000/38000 [==============================] - 4s 103us/step - loss: 0.1559 - acc: 0.9368 - val_loss: 0.1483 - val_acc: 0.9355
    Epoch 4/5
    38000/38000 [==============================] - 4s 102us/step - loss: 0.1362 - acc: 0.9454 - val_loss: 0.1308 - val_acc: 0.9442
    Epoch 5/5
    38000/38000 [==============================] - 4s 114us/step - loss: 0.1228 - acc: 0.9502 - val_loss: 0.1236 - val_acc: 0.9470





    <keras.callbacks.History at 0x7fb0d6031b38>



Looks alright, the training as well as the validation accuracy is climbing from epoch to epoch and slows down a little more after every epoch. It is now your task to find better hyperparameters for our network and training procedure to see how high up you can get the accuracy. For doing that, I suggest taking dl_intro.py, commenting out all the code after the training, adjust your hyperparameters and/or network architectures how you like and running it in the terminal via python dl_intro.py. 

Tipps: 

* Should we train longer?

* Do we need more hidden layers?

* Could we max pool over the entire sequence to get one output per filter?

* Would a second layer be beneficial?

* You could also bias the network architechture with some biological knowledge: How long are transcription factor binding motifs in general and what would be an appropriate filter_width then?

* (Hint: you can do pretty well in 5 - 10 epochs with minor tweaks!)


----

# Evaluation and Prediction

Once you are happy with you network performance or in case you want to jump ahead first and optimize later, we will evaluate our network on the held out test data. Technically, we only optimized on the training data set but we always kept an eye on the validation data loss as well. We are discarding all nets that do well on the training but worse at the validation set (overfitted), therefore we always have an intrinsic bias. The test data set is meant to have never been touched throughout the whole optimization process and we evaluate the perormance of our final model on this set to get an unbiased estimate of its performance.



```python
# Evaluate ---------------------------------------------------------------------
test_file = "./data/pwm_seq_200bp_test_set.txt"
test_labels, test_seqs = read_data(test_file)
score = model.evaluate(test_seqs, test_labels, verbose=0)
print('Test loss:', score[0])
print('Test accuracy:', score[1])
```

    Test loss: 0.123626101106
    Test accuracy: 0.947


Once we are happy with our network we obviously want to employ it as well. Lets say we have a new sequence we want to classify.


```python
# Predictions ------------------------------------------------------------------
# read test sequences again
with open(test_file, "r") as f:
    seqs = []
    labels = []
    for i,l in enumerate(f):
        l = l.rstrip()
        l = l.split("\t")
        seqs.append(l[1])
        labels.append(l[0])

# select a single sequence
single_seq = seqs[0]
single_label = labels[0]

print("Sequence: " + single_seq)

# hot encode        
hotseq = get_hot_coded_seq(single_seq)

# calculate predictions
single_prediction = model.predict(np.expand_dims(hotseq, axis=0))
print("\nClass Prediction \"Probability\":")
print("\tClass 0 = %s" % single_prediction[0][0])
print("\tClass 1 = %s" % single_prediction[0][1])
print("\tClass 2 = %s" % single_prediction[0][2])
print("\tClass 3 = %s" % single_prediction[0][3])

# print the true class
print("\nTrue Class: " + single_label)
```

    Sequence: CATCGTCATATATGTAGTAACCATCCTGTATATATGCCCTCTGGTATATATTGATATATGACCGCCACCTGTTGGCAATATAAAGAGATCTAGTATGCAAGAAGCCCACGACCGAAGTCTCCGCCAGTAGGGGTCAGCCACAAGGGGGCGCTATAGGCGCAATTGCGATACTATATATTACAGCACATGACCACGCGACG
    
    Class Prediction "Probability":
    	Class 0 = 8.92778e-08
    	Class 1 = 0.000697179
    	Class 2 = 0.999115
    	Class 3 = 0.00018805
    
    True Class: 2



```python
# or just run all predictions for 
all_test_predictions = model.predict(test_seqs)
print(all_test_predictions.shape)
print(all_test_predictions[5:8])
```

    (1000, 4)
    [[  7.68844402e-05   4.94505726e-02   7.17186809e-01   2.33285725e-01]
     [  4.68848475e-06   2.78699957e-03   1.99698508e-01   7.97509730e-01]
     [  5.91527879e-01   4.08467114e-01   1.72064711e-06   3.34427318e-06]]


----

# Inspection

Now that we have a reasonably working model, we also want to inspect and see what the net has learned. In applications, we often don't care what the network has learned as long as it performs well and outperforms our competitors. For many research problems however, we are exactly interested in what the network has learned. What features distinguish a cat from a dog or if it comes to decision making (e.g. health care or self driving cars), we obviously want to be able to understand and be able to justify why a certain decision has been chosen and learn how to correct missbehavior.

In genomics we usually want to learn what sequence features distinguish the sequences from one another and map them back to biological properties and factors. The easiest way is to just plot the filter weights. In the first convolutional layer, our filters are just like position weight matrices, multiplying every base at every position with a learned weight and summing the value up (plus a bias and pipe it through the RELU activation function). Unfortunatly, this becomes less straight forward to interpret in deeper layers. There are ways of back engineering and learning the importance of filters in higher layers (e.g. https://github.com/kundajelab/deeplift) but we concern ourself only with the simple first layer here.

We can get the weigths of the filters from the model, save them as .txt files and plot them out. I wrote a wrapper to plot the filter weigths for you in R. Run the code, check the filter_X.txt files and look at the plots and try to interpret them.

* Do any look like transciption factor binding sites you know?
* Do you recognize any sequence features that are not binding motifs?
* Can you simplify the sequences/ motifs from the plot an query them in a transcription factor binding motif database (http://jaspar.genereg.net/)
* What is your best bet: Which sequence motifs did we use for simulating the sequence classes?
* Check the input data. Split them up by class into text files with only the sequences one sequence per line (see example). Query them in standard motif analysis tools (e.g. http://rsat.sb-roscoff.fr/oligo-analysis_form.cgi or http://meme-suite.org/tools/meme). Do these tools find different or similar things?


```python
# Inspect weights --------------------------------------------------------------
model_weights = model.get_weights()
filter_weights = model_weights[0]

# save conv filter weights
for k in range(model_weights[0].shape[2]):
    # save single filter weights
    np.savetxt(("./visualize/filter_%s.txt" % k), filter_weights[:,:,k], delimiter="\t")
```


```bash
%%bash

# Plot them using the supplied R script
Rscript ./helper/plot_sequence_kernel_weights_per_dir.R ./visualize ./visualize plot_weight 10 5
```

     [1] "filter_0.txt" "filter_1.txt" "filter_2.txt" "filter_3.txt" "filter_4.txt"
     [6] "filter_5.txt" "filter_6.txt" "filter_7.txt" "filter_8.txt" "filter_9.txt"
    [1] "Saving Plot filter_0.png"
    [1] "Saving Plot filter_1.png"
    [1] "Saving Plot filter_2.png"
    [1] "Saving Plot filter_3.png"
    [1] "Saving Plot filter_4.png"
    [1] "Saving Plot filter_5.png"
    [1] "Saving Plot filter_6.png"
    [1] "Saving Plot filter_7.png"
    [1] "Saving Plot filter_8.png"
    [1] "Saving Plot filter_9.png"
    [1] "Passing: 10  Skipped: 0"

