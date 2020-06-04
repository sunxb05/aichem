import numpy as np
import os
import random
import shutil
from statistics import mean
from convolutional_neural_network import ConvolutionalNeuralNetwork
import datetime
from logger import Logger


GAMMA = 0.995
MEMORY_SIZE = 100000
BATCH_SIZE = 10
TRAINING_FREQUENCY = 20
# TARGET_NETWORK_UPDATE_FREQUENCY = 40000
TARGET_NETWORK_UPDATE_FREQUENCY = 4000
# MODEL_PERSISTENCE_UPDATE_FREQUENCY = 10000
MODEL_PERSISTENCE_UPDATE_FREQUENCY = 1000
# REPLAY_START_SIZE = 50000

EXPLORATION_MAX = 1.0
EXPLORATION_MIN = 0.1
EXPLORATION_TEST = 0.02
EXPLORATION_DECAY = 0.995



class BaseGameModel:

    def __init__(self, game_name, mode_name, logger_path, observation_space, action_space):
        self.action_space = action_space
        self.observation_space = observation_space
        self.logger = Logger(game_name + " " + mode_name, logger_path)

    def save_run(self, score, step, run):
        self.logger.add_score(score)
        self.logger.add_step(step)
        self.logger.add_run(run)

    def get_move(self, state):
        pass

    def act(self, state):
        pass

    def remember(self, state, action, reward, next_state, done):
        pass

    def step_update(self, total_step):
        pass

    def _get_date(self):
        return str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))


class DDQNGameModel(BaseGameModel):

    def __init__(self, game_name, mode_name, observation_space, action_space, logger_path, model_path):
        BaseGameModel.__init__(self, game_name,
                               mode_name,
                               logger_path,
                               observation_space,
                               action_space)
        self.model_path = model_path
        self.ddqn = ConvolutionalNeuralNetwork(self.observation_space, action_space).model
        if os.path.isfile(self.model_path):
            self.ddqn.load_weights(self.model_path)

    def _save_model(self):
        self.ddqn.save_weights(self.model_path)


class DDQNSolver(DDQNGameModel):

    def __init__(self, game_name, observation_space, action_space):
        testing_model_path = "./output/neural_nets/" + game_name + "/testing/model.h5"
        assert os.path.exists(os.path.dirname(testing_model_path)), "No testing model in: " + str(testing_model_path)
        DDQNGameModel.__init__(self,
                               game_name,
                               "DDQN testing",
                               observation_space,
                               action_space,
                               "./output/logs/" + game_name + "/testing/" + self._get_date() + "/",
                               testing_model_path)

    def act(self, state):
        if np.random.rand() < EXPLORATION_TEST:
            return random.randrange(self.action_space)
        q_values = self.ddqn.predict(state)
        print (np.argmax(q_values[0]))
        return np.argmax(q_values[0])


class DDQNTrainer(DDQNGameModel):

    def __init__(self, game_name, observation_space, action_space):
        DDQNGameModel.__init__(self,
                               game_name,
                               "DDQN training",
                               observation_space,
                               action_space,
                               "./output/logs/" + game_name + "/training/" + self._get_date() + "/",
                               "./output/neural_nets/" + game_name + "/" + self._get_date() + "/model.h5")

        if os.path.exists(os.path.dirname(self.model_path)):
            shutil.rmtree(os.path.dirname(self.model_path), ignore_errors=True)
        os.makedirs(os.path.dirname(self.model_path))

        self.ddqn_target = ConvolutionalNeuralNetwork(self.observation_space, action_space).model
        self._reset_target_network()
        self.exploration_rate = EXPLORATION_MAX
        self.memory = []

    def act(self, state):
        if np.random.rand() < self.exploration_rate:
            return random.randrange(self.action_space)
        q_values = self.ddqn.predict(state)
        return np.argmax(q_values[0])

    def remember(self, state, action, reward, next_state, terminal):
        self.memory.append((state, action, reward, next_state, terminal))
        if len(self.memory) > MEMORY_SIZE:
            self.memory.pop(0)

    def step_update(self, total_step):

        if total_step % TRAINING_FREQUENCY == 0:
            loss, accuracy, average_max_q = self.experience_replay()
            self.logger.add_loss(loss)
            self.logger.add_accuracy(accuracy)
            self.logger.add_q(average_max_q)

        if total_step % MODEL_PERSISTENCE_UPDATE_FREQUENCY == 0:
            self._save_model()

        if total_step % TARGET_NETWORK_UPDATE_FREQUENCY == 0:
            self._reset_target_network()
            print('{{"metric": "total_step", "value": {}}}'.format(total_step))

    def experience_replay(self):
        if len(self.memory) < BATCH_SIZE:
            return
        batch = random.sample(self.memory, BATCH_SIZE)
        for state, action, reward, state_next, terminal in batch:
            q_update = reward
            if not terminal:
                q_update = (reward + GAMMA * np.amax(self.ddqn_target.predict(state_next)[0]))
            q_values = self.ddqn.predict(state)
            q_values[0][action] = q_update

            fit = self.ddqn.fit(state, q_values, batch_size=BATCH_SIZE, verbose=0)
            loss = fit.history["loss"][0]
            accuracy = fit.history["accuracy"][0]
            return loss, accuracy, q_update
        self.exploration_rate *= EXPLORATION_DECAY
        self.exploration_rate = max(EXPLORATION_MIN, self.exploration_rate)

    def _reset_target_network(self):
        self.ddqn_target.set_weights(self.ddqn.get_weights())
