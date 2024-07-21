import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import Sequential
from tensorflow.keras.layers import Dense, Dropout, Activation

class MLP(object):
    def __init__(self, dims):
        self.dims = dims
        self.model = None
        self.input_shape = None
        self.n_classes = None
        self.random_seed = 0 ## for reproducibility

    def init_MLP_model(self, dropout_rate=0.5):
        dense_kernel_init = keras.initializers.TruncatedNormal(mean=0, stddev=0.1, seed=self.random_state) ## same as GEDFN
        model = Sequential()
        model.add(keras.Input(shape=self.input_shape))
        for i in range(len(self.dims)):
            model.add(Dropout(rate=dropout_rate, seed=self.random_state, name="dropout_"+str(i)))
            model.add(Dense(self.dims[i], kernel_initializer=dense_kernel_init, name="dense_"+str(i)))
            model.add(Activation('relu', name="act_"+str(i)))
        model.add(Dense(self.n_classes, kernel_initializer=dense_kernel_init, name="dense_"+str(i+1)))
        self.model = model

    def fit(self, x_train, y_train, batch_size=16, max_epochs=100, 
            sample_weight=None, class_weight=None):
        ## add callback with 5 steps no improvement
        #callback = keras.callbacks.EarlyStopping(monitor='loss', patience=5)
        self.model.fit(x_train, y_train, epochs=max_epochs, batch_size=batch_size, 
                validation_split=0.0, 
                #callbacks=[callback], 
                verbose=2, 
                sample_weight=sample_weight, class_weight=class_weight) # without cross validation

    def compile(self, optimizer='adam'):
        if optimizer == 'adam':
            lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
                            initial_learning_rate=1e-4,
                            decay_steps=1000,
                            decay_rate=0.95, 
                            staircase=True)
            optimizer = tf.keras.optimizers.Adam(learning_rate=lr_schedule)
        self.model.compile(
                loss=keras.losses.CategoricalCrossentropy(from_logits=True),
                metrics=["accuracy"], ## show training accuracy,
                optimizer=optimizer)

    def predict(self, x_test):
        ## with softmax
        x_pred = tf.nn.softmax(self.model.predict(x_test)).numpy()
        return x_pred


class batch_MLP(object):
    def __init__(self, dims):
        self.dims = dims
        self.model = None
        self.input_shape = None
        self.n_batch = None
        self.n_classes = None
        self.random_seed = 0 ## for reproducibility

    def init_batch_MLP_model(self, dropout_rate=0.5):
        dense_kernel_init = keras.initializers.TruncatedNormal(mean=0, stddev=0.1, seed=self.random_state) ## same as GEDFN

        inputs = keras.Input(shape=(self.input_shape))
        batch_inputs = keras.Input(shape=(self.n_batch))
        x = Dropout(rate=dropout_rate, seed=self.random_state)(inputs)
        x = Dense(self.dims[0], kernel_initializer=dense_kernel_init)(x)
        x = Activation('relu')(x)
        batch_x = Dense(self.dims[0], kernel_initializer=dense_kernel_init)(batch_inputs)
        combined = tf.add(x, batch_x)
        for i in range(len(self.dims)-1):
            combined = Dropout(rate=dropout_rate, seed=self.random_state)(combined)
            combined = Dense(self.dims[i+1], kernel_initializer=dense_kernel_init)(combined)
            combined = Activation('relu')(combined)
            batch_combined = Dense(self.dims[i+1], kernel_initializer=dense_kernel_init)(batch_inputs)
            combined = tf.add(combined, batch_combined)
        output = Dense(self.n_classes, kernel_initializer=dense_kernel_init)(combined)
        model = keras.Model(inputs=[inputs, batch_inputs], outputs=output)
        self.model = model

    def fit(self, x_train, b_train, y_train, batch_size=16, max_epochs=100, 
            sample_weight=None, class_weight=None):
        ## add callback with 5 steps no improvement
        #callback = keras.callbacks.EarlyStopping(monitor='loss', patience=5)
        self.model.fit([x_train, b_train], y_train, epochs=max_epochs, batch_size=batch_size, 
                validation_split=0.0, 
                #callbacks=[callback], 
                verbose=2, 
                sample_weight=sample_weight, class_weight=class_weight) # without cross validation

    def compile(self, optimizer='adam'):
        if optimizer == 'adam':
            lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
                            initial_learning_rate=1e-4,
                            decay_steps=1000,
                            decay_rate=0.95, 
                            staircase=True)
            optimizer = tf.keras.optimizers.Adam(learning_rate=lr_schedule)
        self.model.compile(
                loss=keras.losses.CategoricalCrossentropy(from_logits=True),
                metrics=["accuracy"], ## show training accuracy,
                optimizer=optimizer)

    def predict(self, x_test):
        # extract certain layers out
        sub_layers_index = [2, 4, 8, 10]
        included_layers = [self.model.layers[idx] for idx in sub_layers_index]

        inputs = keras.Input(shape=(self.input_shape))
        x = inputs
        for layer in included_layers:
            x = layer(x)
        output = Dense(self.n_classes, kernel_initializer=dense_kernel_init)(x)
        pred_model = keras.Model(inputs=inputs, outputs=output)
        weights_idx = [0, 1, 4, 5, 8, 9]
        original_model_weights = self.model.get_weights()
        pred_model_weights = [original_model_weights[idx] for idx in weights_idx]
        pred_model.set_weights(pred_model_weights)
        ## with softmax
        x_pred = tf.nn.softmax(pred_model.predict(x_test)).numpy()
        return x_pred

