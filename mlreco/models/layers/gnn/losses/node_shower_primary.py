"""Module that defines an EM shower primary identification loss."""

import torch
import numpy as np

from mlreco.models.layers.factories import loss_fn_construct

from mlreco.utils.globals import GROUP_COL, PRGRP_COL
from mlreco.utils.weighting import get_class_weights
from mlreco.utils.gnn.cluster import get_cluster_label_batch
from mlreco.utils.gnn.evaluation import (
        node_assignment_batch, node_assignment_score_batch,
        node_purity_mask_batch)

__all__ = ['NodeShowerPrimaryLoss']


class NodeShowerPrimaryLoss(torch.nn.Module):
    """Loss used to train the EM shower primary identification.

    Takes the two-channel node output of the GNN and optimizes node-wise scores
    such that nodes that initiate a particle cascade are given a high score
    (exclusively relevant for showers for now).

    For use in config:

    ..  code-block:: yaml

        model:
          name: grappa
          modules:
            grappa_loss:
              node_loss:
                name: shower_primary
                <dictionary of arguments to pass to the loss>

    See configuration files prefixed with `grappa_` under the `config`
    directory for detailed examples of working configurations.
    """
    name = 'shower_primary'

    def __init__(self, balance_loss=False, high_purity=False,
                 use_group_pred=False, group_pred_alg='score', loss='ce'):
        """Initialize the EM shower primary identification loss function.

        Parameters
        ----------
        balance_loss : bool, default False
            Whether to weight the loss to account for class imbalance
        high_purity : bool, default False
            Only apply loss to nodes which belong to a sensible group, i.e.
            one with exactly one primary in it (not 0, not > 1)
        use_group_pred : bool, default False
            Use predicted group to check for high purity
        group_pred_alg : str, default 'score'
            Method used to form a predicted group ('threshold' or 'score')
        loss : str, default 'ce'
            Name of the loss function to apply
        """
        # Initialize the parent class
        super().__init__()

        # Initialize basic parameters
        self.balance_loss = balance_loss
        self.high_purity = high_purity
        self.use_group_pred = use_group_pred
        self.group_pred_alg = group_pred_alg

        # Set the loss
        self.loss_fn = loss_fn_construct(loss, functional=True)

    def forward(self, clust_label, clusts, node_pred, edge_index=None,
                edge_pred=None, group_ids=None, **kwargs):
        """Applies the shower primary loss to a batch of data.

        Parameters
        ----------
        clust_label : TensorBatch
            (N, 1 + D + N_f) Tensor of cluster labels for the batch
        clusts : IndexBatch
            (C) Index which maps each cluster to a list of voxel IDs
        node_pred : TensorBatch
            (C, 2) Node prediction logits (binary output)
        edge_index : EdgeIndexBatch, optional
            (2, E) Incidence matrix between clusters
        edge_pred : TensorBatch, optional
            (E, 2) Edge prediction logits (binary output)
        group_ids : TensorBatch, optional
            (C)  Group to which each node belongs to
        **kwargs : dict, optional
            Other labels/outputs of the model which are not relevant here

        Returns
        -------
        loss : torch.Tensor
            Value of the loss
        accuracy : float
            Value of the node-wise classification accuracy
        count : int
            Number of nodes the loss was applied to
        """
        # Fetch the primary and group IDs
        primary_ids = get_cluster_label_batch(
                clust_label, clusts, column=PRGRP_COL)
        if group_ids is None:
            if not self.use_group_pred:
                group_ids = get_cluster_label_batch(
                        clust_label, clusts, column=GROUP_COL)
            else:
                if self.group_pred_alg == 'threshold':
                    edge_pred_np = edge_pred.to_numpy()
                    group_ids = node_assignment_batch(
                            edge_index, edge_pred_np, clusts)

                elif self.group_pred_alg == 'score':
                    edge_pred_np = edge_pred.to_numpy()
                    group_ids = node_assignment_score_batch(
                            edge_index, edge_pred_np, clusts)

                else:
                    raise ValueError("Group prediction algorithm not "
                                     "recognized:", self.group_pred_alg)

        # Create a mask for valid nodes (-1 indicates an invalid primary ID)
        valid_mask = primary_ids.tensor > -1

        # If requested, remove groups that do not contain exactly one primary
        if self.high_purity:
            valid_mask &= node_purity_mask_batch(group_ids, primary_ids)

        # Apply the valid mask and convert the labels to a torch.Tensor
        valid_index = np.where(valid_mask)[0]
        node_assn = primary_ids.to_tensor(
                dtype=torch.long, device=node_pred.device)
        node_assn = node_assn.tensor[valid_index]
        node_pred = node_pred.tensor[valid_index]

        # Compute the loss. Balance classes if requested
        weights = None
        if self.balance_loss:
            weights = get_class_weights(node_assn, num_classes=2)

        loss = self.loss_fn(
                node_pred, node_assn, weight=weights, reduction='sum')
        if len(valid_index):
            loss /= len(valid_index)

        # Compute accuracy of assignment (fraction of correctly assigned nodes)
        acc = 1.
        if len(valid_index):
            acc = float(torch.sum(torch.argmax(node_pred, dim=1) == node_assn))
            acc /= len(valid_index)

        return {
            'accuracy': acc,
            'loss': loss,
            'count': len(valid_index)
        }