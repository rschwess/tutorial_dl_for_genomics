import keras
from keras.models import Sequential
from keras.layers import Dense, Flatten
from keras.layers import Conv1D, MaxPooling1D
import os
import numpy as np

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


# read data --------------------------------------------------------------------
train_file = "./data/pwm_seq_200bp_train_set.txt"
train_labels, train_seqs = read_data(train_file)
test_file = "./data/pwm_seq_200bp_test_set.txt"
test_labels, test_seqs = read_data(test_file)

num_classes = 4

# check data and shapes
# ....

# construct the model ----------------------------------------------------------
model = Sequential()
model.add(Conv1D(10, kernel_size=(10), activation='relu', input_shape=(200, 4), padding='same'))
model.add(MaxPooling1D(pool_size=5))
model.add(Flatten())
model.add(Dense(num_classes, activation='softmax'))
# compile
model.compile(optimizer='adam',
              loss='binary_crossentropy',
              metrics=['accuracy'])

# print model summary
model.summary()

# Train ------------------------------------------------------------------------
model.fit(train_seqs, train_labels,
          batch_size=100,
          epochs=5,
          verbose=1,
          validation_data=(test_seqs, test_labels))


# Evaluate ---------------------------------------------------------------------
valid_file = "./data/pwm_seq_200bp_test_set.txt"
valid_labels, valid_seqs = read_data(test_file)
score = model.evaluate(valid_seqs, valid_labels, verbose=0)
print('Test loss:', score[0])
print('Test accuracy:', score[1])

# Predictions ------------------------------------------------------------------
# calculate predictions
single_prediction = model.predict(np.expand_dims(valid_seqs[1,:], axis=0))
print(single_prediction)
all_valid_predictions = model.predict(valid_seqs)

# Inspect weights --------------------------------------------------------------
model_weights = model.get_weights()
filter_weights = model_weights[0]
# save conv filter weights
for k in range(model_weights[0].shape[2]):
    # save single filter weights
    np.savetxt(("./visualize/filter_%s.txt" % k), filter_weights[:,:,k], delimiter="\t")
# Plot them using the supplied R script
os.system("Rscript ./helper/plot_sequence_kernel_weights_per_dir.R ./visualize ./visualize plot_weight 10 5")
