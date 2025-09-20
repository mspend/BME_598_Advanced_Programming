#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 10:11:24 2025

@author: maurispendlove
"""

# Best visual: bar graph or histogram


import matplotlib.pyplot as plt

# players_dict = {'Adam Chase': 7, 'Ben Doyle': 7, 'Sam Denby': 6, 'Toby Hendy': 2}


# Sample data
players = ['Adam Chase', 'Ben Doyle', 'Sam Denby', 'Toby Hendy']
num_wins = [7,7,6,2]

colors = ['red', 'yellow', 'purple', 'brown']

# Use numeric x values
x_pos = range(len(players))

# Create bar chart
bar_graph = plt.bar(x_pos, num_wins, color=colors)

# Set category labels on the x-axis
plt.xticks(x_pos, players)

# Add title and labels
plt.title('Jet Lag the Game Leaderboard')
plt.xlabel('Players')
plt.ylabel('Number of Wins')
plt.legend()

# Add legend manually
for bar, player in zip(bar_graph, players):
    bar.set_label(player)
plt.legend()

plt.savefig('JetLagLeaderboard.pdf')

# Show the plot
plt.show()


