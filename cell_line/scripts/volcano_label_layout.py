#!/usr/bin/env python3
"""Deterministic, collision-free label centers for a fixed volcano plot panel.

Inputs: JSON containing panel dimensions in physical points, original data-axis
ranges, font path, label identities and unchanged scatter coordinates. Outputs:
JSON containing label centers and leader endpoints in the original data units.
The deterministic search enforces marker visibility and nonoverlapping labels;
plot coordinates, statistics and which genes receive labels never change.

Usage: volcano_label_layout.py input.json output.json
Use '-' for either JSON stream. This moves labels only, never plotted points.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from fontTools.ttLib import TTFont


def font_widths(font_path, strings, font_size):
    """Measure TrueType advances and kerning in points for each complete label."""
    font = TTFont(font_path)
    cmap = font.getBestCmap()
    metrics = font["hmtx"].metrics
    kern = {}
    if "kern" in font:
        for table in font["kern"].kernTables:
            if getattr(table, "coverage", 1) & 1:
                kern.update(getattr(table, "kernTable", {}))
    result = []
    for string in strings:
        if any(ord(char) not in cmap for char in string):
            raise ValueError(f"Font lacks a glyph in {string!r}")
        glyphs = [cmap[ord(char)] for char in string]
        advance = sum(metrics[glyph][0] for glyph in glyphs)
        advance += sum(kern.get(pair, 0) for pair in zip(glyphs, glyphs[1:]))
        result.append(advance * font_size / font["head"].unitsPerEm)
    font.close()
    return np.asarray(result, dtype=float)


def layout(data):
    """Place labels in physical panel space, then convert positions back to data units."""
    # 1. Convert label dimensions and data coordinates into physical points.
    width, height = float(data["panel_width_pt"]), float(data["panel_height_pt"])
    x_range, y_range = np.asarray(data["x_range"], float), np.asarray(data["y_range"], float)
    assert width > 0 and height > 0
    assert x_range.shape == y_range.shape == (2,)
    assert np.isfinite(x_range).all() and np.isfinite(y_range).all()
    assert x_range[1] > x_range[0] and y_range[1] > y_range[0]
    labels = data["labels"]
    identifiers = [r["ensembl_gene_id"] for r in labels]
    strings = [r["display_label"] for r in labels]
    categories = [r["category"] for r in labels]
    assert len(identifiers) == len(set(identifiers)) and all(strings)
    assert all(category in ("Up", "Down") for category in categories)
    n = len(labels)
    font_size = float(data.get("font_size_pt", 8.5))
    margin = 2.0
    clearance = 1.0
    marker_buffer = 1.2
    label_height = font_size * 1.25
    step = float(data.get("grid_step_pt", 2.0))
    assert 0 < step <= 2
    widths = font_widths(data["font_path"], strings, font_size)
    box_widths = widths + 2.0  # 1 pt on each side of the shaped text advance.
    if any(box_widths + 2 * margin > width) or label_height + 2 * margin > height:
        raise ValueError("A label is wider/taller than the panel; increase the panel dimensions.")
    anchors_data = np.asarray([[r["log2FoldChange"], r["minus_log10_padj"]] for r in labels], float).reshape((-1, 2))
    assert np.isfinite(anchors_data).all()
    scale = np.array([width / np.diff(x_range)[0], height / np.diff(y_range)[0]])
    minimum = np.array([x_range[0], y_range[0]])
    anchors = (anchors_data - minimum) * scale

    # A one-point histogram and summed-area table give exact binned density
    # queries without scanning all 17,000+ points for each candidate rectangle.
    points = np.asarray([[p["x"], p["y"]] for p in data.get("points", [])], float).reshape((-1, 2))
    if len(points):
        assert np.isfinite(points).all()
        points = (points - minimum) * scale
    nx, ny = max(1, math.ceil(width)), max(1, math.ceil(height))
    histogram = np.histogram2d(points[:, 1], points[:, 0], bins=[ny, nx], range=[[0, height], [0, width]])[0]
    integral = np.pad(histogram.cumsum(0).cumsum(1), ((1, 0), (1, 0)))

    def density(cx, cy, box_width):
        """Count plotted points inside candidate rectangles using a summed-area table."""
        radius = 0.8
        left = np.clip(np.floor((cx-box_width/2-radius)/width*nx).astype(int), 0, nx)
        right = np.clip(np.ceil((cx+box_width/2+radius)/width*nx).astype(int), 0, nx)
        bottom = np.clip(np.floor((cy-label_height/2-radius)/height*ny).astype(int), 0, ny)
        top = np.clip(np.ceil((cy+label_height/2+radius)/height*ny).astype(int), 0, ny)
        return integral[top, right]-integral[bottom, right]-integral[top, left]+integral[bottom, left]

    # 2. Enumerate deterministic candidates; selected gene markers remain visible.
    candidates, costs, cloud_counts = [], [], []
    endpoints, starts, edge_gaps, ambiguous_flags, draw_flags = [], [], [], [], []
    for i in range(n):
        # Include the eight local label positions around the anchor as well as
        # the complete grid. Their box-edge gap is 2 pt along either axis.
        left, right = margin+box_widths[i]/2, width-margin-box_widths[i]/2
        bottom, top = margin+label_height/2, height-margin-label_height/2
        local_x = anchors[i, 0]+np.array([-box_widths[i]/2-2, 0, box_widths[i]/2+2])
        local_y = anchors[i, 1]+np.array([-label_height/2-2, 0, label_height/2+2])
        xs = np.unique(np.r_[np.arange(left, right+1e-8, step), right, np.clip(local_x, left, right)])
        ys = np.unique(np.r_[np.arange(bottom, top+1e-8, step), top, np.clip(local_y, bottom, top)])
        xx, yy = np.meshgrid(xs, ys)
        xy = np.column_stack([xx.ravel(), yy.ravel()])
        # Keep every selected gene marker visible, including the label's own
        # anchor. Unselected point density remains a soft placement preference.
        anchor_clear = np.ones(len(xy), bool)
        for anchor in anchors:
            anchor_clear &= ((np.abs(xy[:, 0]-anchor[0]) >= box_widths[i]/2+marker_buffer) |
                             (np.abs(xy[:, 1]-anchor[1]) >= label_height/2+marker_buffer))
        xy = xy[anchor_clear]
        if not len(xy):
            raise RuntimeError("Selected markers leave no safe position for a label; enlarge the panel.")
        cloud = density(xy[:, 0], xy[:, 1], box_widths[i])
        center_distance = np.linalg.norm(xy-anchors[i], axis=1)
        box_min = xy-np.array([box_widths[i]/2, label_height/2])
        box_max = xy+np.array([box_widths[i]/2, label_height/2])
        endpoint = np.clip(anchors[i], box_min, box_max)
        gap = np.linalg.norm(endpoint-anchors[i], axis=1)
        assert np.all(gap >= marker_buffer-1e-8)
        ambiguous = np.zeros(len(xy), bool)
        for j, anchor in enumerate(anchors):
            if j != i and categories[j] == categories[i]:
                other_endpoint = np.clip(anchor, box_min, box_max)
                other_gap = np.linalg.norm(other_endpoint-anchor, axis=1)
                ambiguous |= other_gap <= gap+1.5
        draw = (gap > 5.0) | ambiguous
        start_offset = np.minimum(1.5, gap/2)
        start = anchors[i]+(endpoint-anchors[i])*(start_offset/gap)[:, None]
        # Favor adjacent labels, retain necessary association leaders, and make
        # long displacements expensive without turning density into a hard ban.
        cost = (gap+0.05*center_distance+4.0*np.log1p(cloud)+6.0*draw+
                np.maximum(gap-8.0, 0)**2/12.0+2.0*ambiguous)
        candidates.append(xy)
        costs.append(cost)
        cloud_counts.append(cloud)
        endpoints.append(endpoint)
        starts.append(start)
        edge_gaps.append(gap)
        ambiguous_flags.append(ambiguous)
        draw_flags.append(draw)

    def segment_box(x0, y0, x1, y1, left, bottom, right, top):
        """Vectorized closed-segment/rectangle intersection using axis slabs."""
        x0, y0, x1, y1, left, bottom, right, top = np.broadcast_arrays(
            x0, y0, x1, y1, left, bottom, right, top)
        low = np.zeros_like(x0, dtype=float)
        high = np.ones_like(x0, dtype=float)
        for start, end, lower, upper in ((x0, x1, left, right), (y0, y1, bottom, top)):
            delta = end-start
            a = np.divide(lower-start, delta, out=np.full_like(delta, -np.inf, dtype=float), where=delta != 0)
            b = np.divide(upper-start, delta, out=np.full_like(delta, np.inf, dtype=float), where=delta != 0)
            outside = (delta == 0) & ((start < lower) | (start > upper))
            axis_low, axis_high = np.minimum(a, b), np.maximum(a, b)
            axis_low = np.where(outside, np.inf, axis_low)
            axis_high = np.where(outside, -np.inf, axis_high)
            low, high = np.maximum(low, axis_low), np.minimum(high, axis_high)
        return low <= high

    def interaction_score(i, indices, j, other_index):
        """Penalize visible leaders crossing another label or another visible leader."""
        # Penalize only the exact visible segments, from their offset start to
        # the nearest box edge; omitted leaders generate no crossing penalty.
        xy = candidates[i][indices]
        own_start, own_end = starts[i][indices], endpoints[i][indices]
        other = candidates[j][other_index]
        other_start, other_end = starts[j][other_index], endpoints[j][other_index]
        own_draw, other_draw = draw_flags[i][indices], draw_flags[j][other_index]
        candidate_hits = segment_box(own_start[:, 0], own_start[:, 1], own_end[:, 0], own_end[:, 1],
            other[0]-box_widths[j]/2, other[1]-label_height/2,
            other[0]+box_widths[j]/2, other[1]+label_height/2) & own_draw
        existing_hits = segment_box(other_start[0], other_start[1], other_end[0], other_end[1],
            xy[:, 0]-box_widths[i]/2, xy[:, 1]-label_height/2,
            xy[:, 0]+box_widths[i]/2, xy[:, 1]+label_height/2) & other_draw
        def cross(a, b):
            """Return the signed 2D cross product for segment intersection tests."""
            return a[..., 0]*b[..., 1]-a[..., 1]*b[..., 0]
        own_line, other_line = own_end-own_start, other_end-other_start
        first = cross(own_line, other_start-own_start)
        second = cross(own_line, other_end-own_start)
        third = cross(other_line, own_start-other_start)
        fourth = cross(other_line, own_end-other_start)
        crossing = (first*second < -1e-9) & (third*fourth < -1e-9) & own_draw & other_draw
        return 30.0*candidate_hits+30.0*existing_hits+10.0*crossing

    def candidate_score(i, indices, placement):
        """Combine displacement/density cost with interactions against placed labels."""
        score = costs[i][indices].copy()
        for j, candidate in enumerate(placement):
            if j != i and candidate >= 0:
                score += interaction_score(i, indices, j, candidate)
        return score

    def best_candidate(i, placement):
        """Choose a nonoverlapping candidate with deterministic stable tie handling."""
        xy = candidates[i]
        valid = np.ones(len(xy), bool)
        for j, candidate in enumerate(placement):
            if j == i or candidate < 0:
                continue
            other = candidates[j][candidate]
            valid &= ((np.abs(xy[:, 0]-other[0]) >= (box_widths[i]+box_widths[j])/2+clearance-1e-9) |
                      (np.abs(xy[:, 1]-other[1]) >= label_height+clearance-1e-9))
        possible = np.flatnonzero(valid)
        if not len(possible):
            return -1, math.inf
        # Interactions are evaluated on the best 1,200 feasible base-cost
        # positions. All grid positions remain eligible for hard constraints.
        if len(possible) > 1200:
            possible = possible[np.argsort(costs[i][possible], kind="stable")[:1200]]
        score = candidate_score(i, possible, placement)
        best_index = int(np.argmin(score))
        return int(possible[best_index]), float(score[best_index])

    def total_score(placement):
        """Score a complete layout, counting each pairwise interaction only once."""
        total = float(sum(costs[i][placement[i]] for i in range(n)))
        for i in range(n):
            for j in range(i):
                total += float(interaction_score(i, np.array([placement[i]]), j, placement[j])[0])
        return total

    # 3. Try six deterministic placement orders, then improve each by coordinate descent.
    orders = [list(range(n)), sorted(range(n), key=lambda i: (-box_widths[i], i)),
              sorted(range(n), key=lambda i: (-anchors[i, 1], anchors[i, 0], i)),
              sorted(range(n), key=lambda i: (anchors[i, 1], anchors[i, 0], i)),
              sorted(range(n), key=lambda i: (anchors[i, 0], -anchors[i, 1], i)),
              sorted(range(n), key=lambda i: (-anchors[i, 0], -anchors[i, 1], i))]
    best, best_total = None, math.inf
    for order in orders:
        placement = [-1]*n
        for i in order:
            placement[i], _ = best_candidate(i, placement)
            if placement[i] < 0:
                break
        if any(candidate < 0 for candidate in placement):
            continue
        # Deterministic coordinate descent only accepts strict score reductions.
        for _ in range(8):
            changed = False
            improve_order = sorted(range(n), key=lambda i: (-costs[i][placement[i]], i))
            for i in improve_order:
                candidate, score = best_candidate(i, placement)
                current = candidate_score(i, np.array([placement[i]]), placement)[0]
                if candidate >= 0 and score < current-1e-8:
                    placement[i] = candidate
                    changed = True
            if not changed:
                break
        total = total_score(placement)
        if total < best_total-1e-8:
            best, best_total = placement[:], total
    if best is None:
        raise RuntimeError("No collision-free grid layout fits this panel; enlarge the panel or reduce label size.")

    # 4. Require geometric separation before exporting positions in data coordinates.
    centers = np.asarray([candidates[i][best[i]] for i in range(n)]).reshape((-1, 2))
    boxes = np.column_stack([centers[:, 0]-box_widths/2, centers[:, 1]-label_height/2,
                             centers[:, 0]+box_widths/2, centers[:, 1]+label_height/2])
    assert np.all(boxes[:, :2] >= margin-1e-8)
    assert np.all(boxes[:, 2:] <= np.array([width, height])-margin+1e-8)
    for i in range(n):
        assert np.all((np.abs(anchors[:, 0]-centers[i, 0]) >= box_widths[i]/2+marker_buffer-1e-8) |
                      (np.abs(anchors[:, 1]-centers[i, 1]) >= label_height/2+marker_buffer-1e-8))
    gaps = []
    for i in range(n):
        for j in range(i):
            dx = max(0, boxes[i, 0]-boxes[j, 2], boxes[j, 0]-boxes[i, 2])
            dy = max(0, boxes[i, 1]-boxes[j, 3], boxes[j, 1]-boxes[i, 3])
            gaps.append(math.hypot(dx, dy))
            assert max(dx, dy) >= clearance-1e-8

    result = []
    for i, original in enumerate(labels):
        center, anchor = centers[i], anchors[i]
        endpoint, start = endpoints[i][best[i]], starts[i][best[i]]
        gap = float(edge_gaps[i][best[i]])
        ambiguous, draw = bool(ambiguous_flags[i][best[i]]), bool(draw_flags[i][best[i]])
        assert draw == (gap > 5.0 or ambiguous)
        label_data = center/scale+minimum
        leader_data = endpoint/scale+minimum
        leader_start_data = start/scale+minimum
        result.append({"id": identifiers[i], "ensembl_gene_id": identifiers[i],
            "display_label": strings[i], "category": categories[i],
            "label_x": float(label_data[0]), "label_y": float(label_data[1]),
            "leader_x": float(leader_data[0]), "leader_y": float(leader_data[1]), "draw_leader": draw,
            "leader_start_x": float(leader_start_data[0]), "leader_start_y": float(leader_start_data[1]),
            "leader_length_pt": gap, "leader_ambiguous": ambiguous,
            "drawn_segment_length_pt": float(np.linalg.norm(endpoint-start)) if draw else 0.0,
            "center_pt": center.tolist(), "bbox_pt": boxes[i].tolist(), "anchor_pt": anchor.tolist(),
            "text_width_pt": float(widths[i]), "reserved_width_pt": float(box_widths[i]),
            "reserved_height_pt": label_height, "nearby_point_count": int(cloud_counts[i][best[i]])})
    return {"labels": result}


def main():
    """Read JSON from a path or stdin and write positioned labels to a path or stdout."""
    if len(sys.argv) != 3:
        raise SystemExit("Usage: volcano_label_layout.py input.json output.json ('-' for standard streams)")
    text = sys.stdin.read() if sys.argv[1] == "-" else Path(sys.argv[1]).read_text()
    result = layout(json.loads(text))
    encoded = json.dumps(result, indent=2, allow_nan=False)+"\n"
    if sys.argv[2] == "-":
        sys.stdout.write(encoded)
    else:
        Path(sys.argv[2]).write_text(encoded)


if __name__ == "__main__":
    main()
