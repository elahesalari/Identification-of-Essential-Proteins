import numpy as np
import pandas as pd
import networkx as nx
from nltk import flatten
import math
from sklearn.metrics import confusion_matrix


class EssentialProteins:
    def __init__(self, gama, max_iteration):
        self.gama = gama
        self.max_iteration = max_iteration

    def read_data(self):
        self.spellman = pd.read_csv('Datasets/Spellman.csv')
        with open('Datasets/Gavin_PPI.txt') as f:
            self.gavin_ppi = f.read()
        with open('Datasets/Essential Genes.txt') as f:
            self.essential = f.read()

    def graph(self):
        max_data = np.max(np.array(self.spellman.iloc[:, 1:]))
        min_data = np.min(np.array(self.spellman.iloc[:, 1:]))

        # calculate ECC
        genes_2 = [[s.strip() for s in line.split('\t') if s] for line in self.gavin_ppi.split('\n') if line]
        flatten_genes = flatten(genes_2)
        unique_genes = np.unique(flatten_genes)

        # Create Graph
        G = nx.DiGraph()
        G.add_nodes_from(unique_genes)
        G.add_edges_from(genes_2)

        # Calculate ECC
        ECC = {}
        for edge in G.edges():
            gn_U = list(G.neighbors(edge[0]))
            gn_V = list(G.neighbors(edge[1]))
            intersection = list(set(gn_U).intersection(set(gn_V)))
            count_intersection = len(intersection)
            degree_U = G.degree(edge[0])
            degree_V = G.degree(edge[1])
            weight = (count_intersection + 1) / min(degree_U, degree_V)
            ECC[(edge[0], edge[1])] = weight

        # Add weight to each edges of graph
        nx.set_edge_attributes(G, values=ECC, name='weight')

        # Calculate PCC
        PCC = {}
        opposite_edges = []
        for edge in G.edges():
            opposite_edges.append([edge[1], edge[0]])
            g_u = self.spellman.loc[self.spellman['time'] == edge[0]]
            g_v = self.spellman.loc[self.spellman['time'] == edge[1]]

            if g_u.empty:
                PCC[(edge[1], edge[0])] = 0
                continue
            if g_v.empty:
                PCC[(edge[1], edge[0])] = 0
                continue
            g_u = g_u.iloc[:, 1:].values.tolist()
            g_v = g_v.iloc[:, 1:].values.tolist()
            g_u = flatten(g_u)
            g_v = flatten(g_v)

            # Normalize data
            g_u = (g_u - min_data) / (max_data - min_data)
            g_v = (g_v - min_data) / (max_data - min_data)

            g_hat_u = np.mean(g_u)
            g_hat_v = np.mean(g_v)

            numerator = np.sum((g_v - g_hat_v) * (g_u - g_hat_u))
            denominator = np.sqrt(np.sum((g_v - g_hat_v) ** 2) * np.sum((g_u - g_hat_u) ** 2))
            weight = numerator / denominator

            PCC[(edge[1], edge[0])] = weight

        G.add_edges_from(opposite_edges)
        nx.set_edge_attributes(G, values=PCC, name='weight')

        return G

    def HSEP(self, G):

        authority = {}
        hub = {}
        for node in G.nodes():
            hub[node] = 1
            authority[node] = 1

        for m in range(1, self.max_iteration):
            print(m)
            a_m_1 = authority.copy()
            h_m_1 = hub.copy()

            for node in G.nodes():

                list_B = [(ed[0], ed[2]['weight']) for ed in G.edges(data=True) if ed[1] == node]
                list_F = [(ed[1], ed[2]['weight']) for ed in G.edges(data=True) if ed[0] == node]

                a_v = 0
                for u in list_B:
                    a_v += hub[u[0]] * u[1]
                authority[node] = a_v
                h_v = 0
                for u in list_F:
                    h_v += authority[u[0]] * u[1]
                hub[node] = h_v

            # Normalization Authority and Hub
            max_a = np.max([v for k, v in authority.items()])
            max_h = np.max([v for k, v in hub.items()])

            for k, v in authority.items():
                authority[k] = v / max_a
            for k, v in hub.items():
                hub[k] = v / max_h

            # Final condition
            for node in G.nodes():
                a_m = authority[node]
                h_m = hub[node]

                if (abs(a_m - a_m_1[node]) + abs(h_m - h_m_1[node])) < self.gama:
                    return authority, hub

    def calc_HSEP_score(self, G, authority, hub):
        alpha = np.linspace(0, 1, 11)
        score = np.zeros((len(alpha), len(G.nodes())))
        k = len(alpha)
        N = [0.01, 0.05, 0.10, 0.15, 0.20, 0.25]
        top_n = []
        for i in N:
            top_n.append(int(len(G.nodes()) * i))
        print(top_n)
        for n in range(len(top_n)):
            top_ranked_in_each_alpha = []
            for i in range(len(alpha)):
                nodes_name = []
                for j, node in enumerate(G.nodes()):
                    nodes_name.append([j, node])
                    score[i, j] = (alpha[i] * authority[node]) + ((1 - alpha[i]) * hub[node])

                # Sort Values and give n top rank from this
                sort_scores = np.sort(score[i])
                rank = sort_scores[-top_n[n]:]

                # Find n top ranked gene names
                indexes = [score[i].tolist().index(r) for r in rank]
                top_ranked_genes = []
                for idx in indexes:
                    for pair in nodes_name:
                        if idx == pair[0]:
                            top_ranked_genes.append(pair[1])
                            break
                top_ranked_in_each_alpha.append(top_ranked_genes)

            self.ensemble_score(top_ranked_in_each_alpha, k, N[n])

    def ensemble_score(self, top_ranked, k, n):
        threshold = np.floor(k / 2)
        EM = set()
        essential_protein = []
        nonessential_protein = []
        for node in G.nodes():
            count = sum(X.count(node) for X in top_ranked)
            EM.add((node, count))
            if count > threshold:
                essential_protein.append(node)
            if count <= threshold:
                nonessential_protein.append(node)

        essential_protein = np.array(list(zip(essential_protein, [1] * len(essential_protein))))
        nonessential_protein = np.array(list(zip(nonessential_protein, [0] * len(nonessential_protein))))
        predict_genes = np.concatenate((essential_protein, nonessential_protein))

        self.evaluate_performance(predict_genes, n)

    def evaluate_performance(self, predict_genes, n):
        print(f'------------*** Top {n * 100}% ***------------')
        dataset_essential = [line.strip() for line in self.essential.split('\n')]

        true_genes = np.zeros((len(G.nodes), 2), dtype=object)
        for i, node in enumerate(G.nodes()):
            true_genes[i, 0] = node
            if node in dataset_essential:
                true_genes[i, 1] = 1

        y_true = true_genes[:, 1].tolist()
        y_predict = list(map(int, predict_genes[:, 1]))

        TN, FP, FN, TP = confusion_matrix(y_true, y_predict).ravel()

        SN = TP / (TP + FN)
        SP = TN / (TN + FP)
        PPV = TP / (TP + FP)
        NPV = TN / (TN + FN)
        F_measure = (2 * SN * PPV) / (SN + PPV)
        ACC = (TP + TN) / (TP + TN + FP + FN)
        print('SN:', SN)
        print('SP:', SP)
        print('PPV:', PPV)
        print('NPV:', NPV)
        print('F-measure:', F_measure)
        print('Accuracy:', ACC)


if __name__ == '__main__':
    maxiter = 10
    gama = 0.002
    protein = EssentialProteins(gama, maxiter)
    protein.read_data()
    G = protein.graph()
    auth, hub = protein.HSEP(G)
    protein.calc_HSEP_score(G, auth, hub)
