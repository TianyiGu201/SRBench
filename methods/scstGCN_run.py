# SPDX-License-Identifier: GPL-3.0-only
"""Run scstGCN in its own installed environment and upstream repository.

Adapted from the scstGCN workflow (GPL-3.0); see THIRD_PARTY_NOTICES.md.
"""
import argparse
import os
from pathlib import Path
import pickle
import sys

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work_dir', type=Path)
    parser.add_argument('--method-dir', type=Path, default=Path.cwd())
    parser.add_argument('--uni-checkpoint', type=Path, default=os.environ.get('UNI_CHECKPOINT'))
    args = parser.parse_args(argv)
    if args.uni_checkpoint is None or not args.uni_checkpoint.expanduser().is_file():
        parser.error('Supply --uni-checkpoint or set UNI_CHECKPOINT to a local checkpoint file')
    os.environ['UNI_CHECKPOINT'] = str(args.uni_checkpoint.expanduser().resolve())
    from install_patches import verify_patches
    verify_patches('scstGCN', args.method_dir)
    sys.path.insert(0, str(args.method_dir.resolve()))
    import align_image_pos
    from histological_feature_map import main as his_ex
    from pos_feature_map import main as pos_ex
    from RGB_feature_map import main as rgb_ex
    import get_mask
    from utils import get_HVG_genes
    from impute import main as train_model
    prefix = os.path.join(str(args.work_dir.resolve()), '')
    align_image_pos.main(prefix)

    login_token = os.environ.get("HF_TOKEN", "")
    img_emb = his_ex(prefix, login_token)
    rgb_emb = rgb_ex(prefix)
    pos_emb = pos_ex(prefix, img_emb)
    embs = dict(his=img_emb, rgb=rgb_emb, pos=pos_emb)
    with open(prefix + 'Multimodal_feature_map.pickle', 'wb') as file:
        pickle.dump(embs, file)


    get_mask.main(prefix)
    if not os.path.isfile(os.path.join(prefix, "gene-names.txt")):
        get_HVG_genes(prefix, n_top=1000)

    train_model(prefix, epoch=500, device='cuda', load_saved=False)

if __name__ == '__main__':
    main()
