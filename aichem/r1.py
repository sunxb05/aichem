import random
from aichem_env import ChemEnv
import numpy as np
from collections import deque
from scores.score_logger import ScoreLogger
from convolutional_neural_network import ConvolutionalNeuralNetwork


GAMMA = 0.95
LEARNING_RATE = 0.001

MEMORY_SIZE = 1000000
BATCH_SIZE = 5

EXPLORATION_MAX = 1.0
EXPLORATION_MIN = 0.01
EXPLORATION_DECAY = 0.995


class DQNActor:

    def __init__(self, observation_space, action_space):
        self.exploration_rate = EXPLORATION_MAX
        self.action_space = action_space
        self.memory = deque(maxlen=MEMORY_SIZE)
        self.model = ConvolutionalNeuralNetwork(observation_space, action_space).model

    def remember(self, state, action, reward, next_state, terminal):
        self.memory.append((state, action, reward, next_state, terminal))

    def act(self, state):
        #  select either value
        #  np.random.rand()     random number between 0–1
        #  random.randrange()   random number (always 0 in this situation)
        if np.random.rand() < self.exploration_rate:
            return random.randrange(self.action_space)
        #  q_values = [[7.9756827 7.7895737]]
        #  np.argmax(q_values[0]) = 0.
        q_values = self.model.predict(state)
        return np.argmax(q_values[0])

    def experience_replay(self):
        if len(self.memory) < BATCH_SIZE:
            return
        batch = random.sample(self.memory, BATCH_SIZE)
        for state, action, reward, state_next, terminal in batch:
            # reward = -1
            q_update = reward
            # reward = 1
            # np.amax  The maximum value along a given axis.
            if not terminal:
                q_update = (reward + GAMMA * np.amax(self.model.predict(state_next)[0]))
            q_values = self.model.predict(state)
            q_values[0][action] = q_update
            self.model.fit(state, q_values, verbose=0)
        self.exploration_rate *= EXPLORATION_DECAY
        self.exploration_rate = max(EXPLORATION_MIN, self.exploration_rate)

def rl():
    env = ChemEnv()
    score_logger = ScoreLogger('aichem')
    # observation_space = env.observation_space.shape[0]
    # action_space = env.action_space.n
    observation_space = 4
    action_space      = 2
    actor = DQNActor(observation_space, action_space)
    run = 0
    while True:
        run += 1
        state = env.reset()
        state = np.reshape(state, [1, observation_space])
        step = 0
        while True:
            step += 1
            action = actor.act(state)
            step_name = str(run)+"_"+str(step)
            state_next, reward, terminal, info = env.step(action, step_name)
            reward = reward if not terminal else -reward
            state_next = np.reshape(state_next, [1, observation_space])
            actor.remember(state, action, reward, state_next, terminal)
            state = state_next
            if terminal:
                print ("Run: " + str(run) + ", exploration: " + str(actor.exploration_rate) + ", score: " + str(step))
                score_logger.add_score(step, run)
                break
            actor.experience_replay()

if __name__ == "__main__":
    rl()
