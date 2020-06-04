from keras.models import Sequential
from keras.layers import Dense
from keras.optimizers import Adam
from collections import deque
LEARNING_RATE = 0.001

MEMORY_SIZE = 1000000
EXPLORATION_MAX = 1.0
EXPLORATION_MIN = 0.01

class ConvolutionalNeuralNetwork:

    def __init__(self, observation_space, action_space):
        self.exploration_rate = EXPLORATION_MAX

        self.action_space = action_space
        self.memory = deque(maxlen=MEMORY_SIZE)

        self.model = Sequential()
        self.model.add(Dense(24, input_shape=(observation_space,), activation="relu"))
        self.model.add(Dense(24, activation="relu"))
        self.model.add(Dense(self.action_space, activation="linear"))
        self.model.compile(loss="mse",
                           optimizer=Adam(lr=LEARNING_RATE),
                           metrics=["accuracy"])
        self.model.summary()
