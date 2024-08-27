#!/usr/bin/env python3

import sys
import argparse
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import plotly.express as px
import plotly

def main():
    parser = argparse.ArgumentParser(description='Recluster all data points using PyTorch.')
    parser.add_argument('--df', required=True, type=str, help='CSV file with the clustered data.')
    parser.add_argument('--confidence', type=float, default=0.9, help='Confidence threshold for assigning clusters.')
    parser.add_argument('--write', action='store_true', help='Flag to write the updated DataFrame to the CSV file.')
    parser.add_argument('--target_accuracy', type=float, required=True, help='Target accuracy for model optimization.')
    parser.add_argument('--test_size', type=float, default=0.2, help='Proportion of the data to use as test set.')
    args = parser.parse_args()
    
    # Load the data
    proteins = pd.read_csv(args.df)
    proteins.set_index('Protein Name', inplace=True)
    # Prepare the data for PyTorch
    classified_data, unclassified_data = prepare_data(proteins)
    
    if classified_data.empty:
        print("No classified data points to train on.")
        return
    
    # Split the data into training and test sets
    train_data, test_data = train_test_split(classified_data, test_size=args.test_size, random_state=42)
    
    best_model, best_scaler, best_accuracy, best_epoch, best_batch_size, best_hidden_size = optimize_model(
        train_data, test_data, args.target_accuracy
    )
    
    # Predict all data points
    predictions = predict_all(best_model, best_scaler, proteins, args.confidence)
    
    # Count how many unclassified points have been classified and how many classified points changed labels
    original_clusters = proteins['Cluster'].astype(int).copy()
    proteins['Cluster'] = predictions
    
    # Count the changes
    start_unclassified_count = np.sum(original_clusters == -1)
    end_unclassified_count = np.sum(proteins['Cluster'] == -1)
    reclassified_count = np.sum((original_clusters != proteins['Cluster'].astype(int)) & (original_clusters != -1))
    
    # Convert cluster labels back to string
    proteins['Cluster'] = proteins['Cluster'].astype(str)

    # Output counts
    print(f"Number of unclassified points start: {start_unclassified_count}")
    print(f"Number of unclassified points end: {end_unclassified_count}")
    print(f"Number of previously classified points that changed labels: {reclassified_count}")

    # Output the final parameters and accuracy
    print(f"Best model parameters - Epochs: {best_epoch}, Batch size: {best_batch_size}, Hidden size: {best_hidden_size}")
    print(f"Best model accuracy on test set: {best_accuracy:.4f}")

    # Write to CSV if --write flag is specified
    if args.write:
        proteins.to_csv(args.df, index=False)
        print("Updated DataFrame written to the CSV file.")
    else:
        print("Run completed. Use --write flag to save changes.")
    
    # Plot the updated data
    fig = plot(proteins)
    plotly.offline.plot(fig, filename='../outputs/dbscan_reclustered_plot.html', auto_open=True)

def prepare_data(proteins):
    # Convert cluster labels to integers, filtering out '-1'
    proteins['Cluster'] = proteins['Cluster'].astype(int)
    classified_data = proteins[proteins['Cluster'] != -1]
    unclassified_data = proteins[proteins['Cluster'] == -1]
    
    return classified_data, unclassified_data

class ClusterNet(nn.Module):
    def __init__(self, input_size, hidden_size, num_classes):
        super(ClusterNet, self).__init__()
        self.fc1 = nn.Linear(input_size, hidden_size)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(hidden_size, hidden_size // 2)
        self.fc3 = nn.Linear(hidden_size // 2, num_classes)
    
    def forward(self, x):
        out = self.fc1(x)
        out = self.relu(out)
        out = self.fc2(out)
        out = self.relu(out)
        out = self.fc3(out)
        return out

def train_model(data, epochs, batch_size, hidden_size, lr):
    features = data[['New Length', 'New pI', 'New Instability', 'New Gravy']].values
    labels = data['Cluster'].values
    
    # Standardize the data
    scaler = StandardScaler()
    features = scaler.fit_transform(features)
    
    # Convert to PyTorch tensors
    features = torch.tensor(features, dtype=torch.float32)
    labels = torch.tensor(labels, dtype=torch.long)
    
    # Build the model
    input_size = features.shape[1]
    num_classes = len(np.unique(labels))  # number of classes excluding -1
    model = ClusterNet(input_size, hidden_size, num_classes)
    
    # Define loss and optimizer
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(model.parameters(), lr=lr)
    
    # Train the model
    for epoch in range(epochs):
        model.train()
        optimizer.zero_grad()
        outputs = model(features)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()
        
        # if (epoch + 1) % 10 == 0:
        #     print(f'Epoch [{epoch + 1}/{epochs}], Loss: {loss.item():.4f}')
    
    return model, scaler

def evaluate_model(model, scaler, test_data):
    features = test_data[['New Length', 'New pI', 'New Instability', 'New Gravy']].values
    labels = test_data['Cluster'].values
    
    # Standardize the data
    features = scaler.transform(features)
    
    # Convert to PyTorch tensors
    features = torch.tensor(features, dtype=torch.float32)
    labels = torch.tensor(labels, dtype=torch.long)
    
    # Predict
    model.eval()
    with torch.no_grad():
        outputs = model(features)
        _, predictions = torch.max(outputs, 1)
    
    # Calculate accuracy
    accuracy = (predictions == labels).sum().item() / labels.size(0)
    
    return accuracy

def optimize_model(train_data, test_data, target_accuracy):
    best_model = None
    best_scaler = None
    best_accuracy = 0
    best_epoch = 0
    best_batch_size = 0
    best_hidden_size = 0
    closest_accuracy = float('inf')
    epochs = [1,5,10,50,100,500]
    batch_sizes = [4,8,16,32,64]
    hidden_sizes = [4,8,16,32,64]
    lrs = [0.00001,0.0001,0.001,0.01,0.1]
    total = len(epochs)*len(batch_sizes)*len(hidden_sizes)*len(lrs)
    count = 0
    for epoch in epochs:
        for batch_size in batch_sizes:
            for hidden_size in hidden_sizes:
                for lr in lrs:
                    count += 1
                    model, scaler = train_model(train_data, epoch, batch_size, hidden_size, lr)
                    accuracy = evaluate_model(model, scaler, test_data)
                    
                    if abs(accuracy - target_accuracy) < closest_accuracy:
                        closest_accuracy = abs(accuracy - target_accuracy)
                        best_model = model
                        best_scaler = scaler
                        best_accuracy = accuracy
                        best_epoch = epoch
                        best_batch_size = batch_size
                        best_hidden_size = hidden_size
                    sys.stdout.flush()
                    percent = round(count*100/total, 0)
                    sys.stdout.write("\rPecent done: {}%, accuracy: {}".format(str(percent), str(best_accuracy)))
               
    return best_model, best_scaler, best_accuracy, best_epoch, best_batch_size, best_hidden_size

def predict_all(model, scaler, data, confidence_threshold):
    features = data[['New Length', 'New pI', 'New Instability', 'New Gravy']].values
    features = scaler.transform(features)
    
    # Convert to PyTorch tensor
    features = torch.tensor(features, dtype=torch.float32)
    
    # Predict
    model.eval()
    with torch.no_grad():
        outputs = model(features)
        confidences, predictions = torch.max(nn.functional.softmax(outputs, dim=1), 1)
    
    # Apply confidence threshold
    predictions = predictions.numpy()
    confidences = confidences.numpy()
    predictions = np.where(confidences >= confidence_threshold, predictions, -1)
    
    return predictions

def color_palette():
    global tab_20, tab_20_pastel
    old_tab_20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 179, 71),  
         (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
         (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
         (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
         (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
    tab_20_pastel = []
    tab_20 = []
    counter = 0
    for color in old_tab_20:
        if counter % 2 == 1:
            tab_20_pastel.append(f'rgb{color}')
        else:
            tab_20.append(f'rgb{color}')
        counter += 1

def plot(proteins):
    color_palette()
    
    # Sort the DataFrame by 'Cluster' in numerical order
    proteins['Cluster'] = proteins['Cluster'].astype(str)
    proteins.sort_values('Cluster', axis=0, ascending=True, inplace=True, kind='quicksort')
    
    # Plot
    fig = px.scatter_3d(proteins, x='Length', y='pI', z='Gravy', color='Cluster', color_discrete_sequence=tab_20_pastel, 
                        hover_name=proteins.index, hover_data=['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum'])
    fig.update_traces(marker=dict(size=3, opacity=0.6))
    fig.update_layout(title={'text': 'PyTorch reclustering of unlabeled points'}, legend={'itemsizing': 'constant'}) 
    fig.update_scenes(camera_eye_x=-1.5)
    fig.update_scenes(camera_eye_y=-1.75)
    fig.update_scenes(camera_eye_z=0.15)
    fig.update_layout(width=1000, height=750)
    return fig

if __name__ == "__main__":
    main()
