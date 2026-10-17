# test_tracking.py
import numpy as np
from scipy import ndimage as ndi


import os.path
import sys

dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, "../libpy/")
sys.path.append(filename)
from breaking import btrack


S8 = ndi.generate_binary_structure(2, 2)
NY, NX = 64, 64
RESULTS = []  # (title, frames, ids), filled by run()


def blob_frame(centres, radius=3):
    """Disc(s) on a periodic grid, centres in cells (cy, cx)."""
    Y, X = np.mgrid[0:NY, 0:NX]
    m = np.zeros((NY, NX), bool)
    for cy, cx in centres:
        dy = (Y - cy + NY / 2) % NY - NY / 2
        dx = (X - cx + NX / 2) % NX - NX / 2
        m |= dy**2 + dx**2 <= radius**2
    return m


def advected(starts, v, nt):
    return [
        blob_frame([(y + v[0] * t, x + v[1] * t) for y, x in starts]) for t in range(nt)
    ]


def unique_ids(ids):
    return np.unique(np.concatenate([np.unique(a[a > 0]) for a in ids]))


def run(title, starts, v, nt, r_max):
    """Build an advected-blob case, track it, and store it for plotting."""
    frames = advected(starts, v, nt)
    ids = btrack(frames, 1.0, 1.0, S8, r_max=r_max)
    RESULTS.append((title, frames, ids))
    return ids


def circ_mean(v, L):
    """Mean of positions on a periodic axis of length L."""
    a = 2 * np.pi * np.asarray(v) / L
    return (
        (np.arctan2(np.sin(a).mean(), np.cos(a).mean()) % (2 * np.pi)) * L / (2 * np.pi)
    )


def centroids(ids):
    """{track id: array of (t, cy, cx)} using circular means."""
    out = {}
    for t, a in enumerate(ids):
        for i in np.unique(a[a > 0]):
            ys, xs = np.nonzero(a == i)
            out.setdefault(i, []).append((t, circ_mean(ys, NY), circ_mean(xs, NX)))
    return {i: np.array(v) for i, v in out.items()}


def plot_results(n_snap=5):
    import matplotlib.pyplot as plt

    for title, frames, ids in RESULTS:
        nt = len(ids)
        times = np.unique(np.linspace(0, nt - 1, n_snap).astype(int))
        fig, ax = plt.subplots(
            1,
            len(times) + 1,
            figsize=(3 * (len(times) + 1), 3.2),
            constrained_layout=True,
        )
        for a, t in zip(ax, times):
            img = np.ma.masked_equal(ids[t], 0)
            a.imshow(
                np.where(img.mask, np.nan, (img.data - 1) % 10),
                origin="lower",
                cmap="tab10",
                vmin=0,
                vmax=9,
                interpolation="nearest",
                extent=(0, NX, 0, NY),
            )
            for i, c in centroids([ids[t]]).items():
                a.text(c[0, 2], c[0, 1], str(i), ha="center", va="center", fontsize=9)
            a.set_title(f"t={t}")
            a.set_xlabel("x (cells)")
        ax[0].set_ylabel("y (cells)")

        # trajectories coloured by time, broken at periodic seams
        from matplotlib.collections import LineCollection
        from matplotlib.cm import ScalarMappable

        a = ax[-1]
        norm = plt.Normalize(0, nt - 1)
        for i, c in centroids(ids).items():
            t, cy, cx = c.T
            pts = np.column_stack([cx, cy])
            segs = np.stack([pts[:-1], pts[1:]], axis=1)
            jump = (np.abs(np.diff(cx)) > NX / 2) | (np.abs(np.diff(cy)) > NY / 2)
            segs = segs[~jump]  # drop segments crossing a seam
            tmid = (0.5 * (t[:-1] + t[1:]))[~jump]
            lc = LineCollection(segs, cmap="viridis", norm=norm, linewidths=2)
            lc.set_array(tmid)
            a.add_collection(lc)
            a.scatter(cx, cy, c=t, cmap="viridis", norm=norm, s=12, zorder=3)
            a.plot(cx[0], cy[0], "k*", ms=9, zorder=4)  # star = start
            a.text(cx[0], cy[0], f" {i}", fontsize=8, va="bottom")  # track id
        a.set(xlim=(0, NX), ylim=(0, NY), title="trajectories", xlabel="x (cells)")
        a.set_aspect(1)
        fig.colorbar(
            ScalarMappable(norm=norm, cmap="viridis"),
            ax=a,
            label="time step",
            shrink=0.8,
        )

        # if len(centroids(ids)) <= 8:
        #     a.legend(fontsize=7, loc="upper right")
        fig.suptitle(title)
    plt.show()


def run_frames(title, frames, r_max):
    events = []
    ids = btrack(frames, 1.0, 1.0, S8, r_max=r_max, events=events)
    RESULTS.append((title, frames, ids))
    return ids, events


def split_frames(y0=32, nt=12, ts=4, vx=2, sep=2):
    """One disc moving in +x. After frame ts it divides in two along y,
    each half moving away by `sep` cells per frame."""
    frames = []
    for t in range(nt):
        x = 10 + vx * t
        if t <= ts:
            centres = [(y0, x)]
        else:
            off = sep * (t - ts)
            centres = [(y0 - off, x), (y0 + off, x)]
        frames.append(blob_frame(centres))
    return frames


def n_ids(a):
    return len(np.unique(a[a > 0]))


def check_split(y0):
    ids, ev = run_frames(f"split, y0={y0}", split_frames(y0), r_max=8.0)
    assert n_ids(ids[5]) == 1 and n_ids(ids[6]) == 2, "frames do not split as designed"
    splits = [e for e in ev if e[0] == "split"]
    assert len(splits) == 1, ev
    _, t, child, parent = splits[0]
    assert (t, parent, child) == (6, 1, 2), splits
    assert set(np.unique(ids[-1][ids[-1] > 0])) == {1, 2}  # parent id survives
    assert [e for e in ev if e[0] == "birth"] == [("birth", 0, 1)]
    assert not [e for e in ev if e[0] in ("merge", "death")], ev


def test_split_records_parentage():
    check_split(32)


def test_split_across_y_seam():
    check_split(2)  # children wrap around y = 0


def test_merge_ends_one_id_and_records_it():
    ids, ev = run_frames("merge (reversed split)", split_frames()[::-1], r_max=8.0)
    assert n_ids(ids[5]) == 2 and n_ids(ids[6]) == 1, "frames do not merge as designed"
    merges = [e for e in ev if e[0] == "merge"]
    assert len(merges) == 1, ev
    _, t, ended, into = merges[0]
    assert t == 6 and {ended, into} == {1, 2}, merges
    assert list(np.unique(ids[6][ids[6] > 0])) == [into]
    assert all(ended not in a for a in ids[6:])  # ended id never returns
    assert not [e for e in ev if e[0] in ("split", "death")], ev
    assert [e[0] for e in ev if e[1] > 0] == ["merge"]  # nothing else after t=0


def check_single_blob(v):
    ids = run(f"single blob, v={v}", [(32, 10)], v, 40, r_max=8.0)
    assert list(unique_ids(ids)) == [1], f"ids found: {unique_ids(ids)}"
    assert all(len(np.unique(a[a > 0])) == 1 for a in ids)


def test_two_blobs_keep_separate_ids():
    ids = run("two blobs, v=(0, 5)", [(15, 10), (45, 40)], (0, 5), 30, r_max=8.0)
    for a in ids:
        assert len(np.unique(a[a > 0])) == 2
    assert ids[0][15, 10] == ids[-1][15, (10 + 5 * 29) % NX]
    assert ids[0][45, 40] == ids[-1][45, (40 + 5 * 29) % NX]
    assert ids[0][15, 10] != ids[0][45, 40]


def test_r_max_too_small_breaks_the_track():
    ids = run(
        "r_max=2, speed 12 (track should break)", [(32, 10)], (0, 12), 10, r_max=2.0
    )
    assert len(unique_ids(ids)) > 1


def test_large_r_max_bridges_the_same_jump():
    ids = run(
        "r_max=8, speed 12 (track should survive)", [(32, 10)], (0, 12), 10, r_max=8.0
    )
    assert list(unique_ids(ids)) == [1]


def test_single_blob_plus_x():
    check_single_blob((0, 5))


def test_single_blob_minus_x():
    check_single_blob((0, -5))


def test_single_blob_plus_y():
    check_single_blob((4, 0))


def test_single_blob_diagonal():
    check_single_blob((3, 4))


if __name__ == "__main__":
    tests = [
        (name, f)
        for name, f in sorted(globals().items())
        if name.startswith("test_") and callable(f)
    ]
    failed = 0
    for name, f in tests:
        try:
            f()
            print(f"PASS  {name}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL  {name}: {e}")
        except Exception as e:  # unexpected error, not a failed assertion
            failed += 1
            print(f"ERROR {name}: {type(e).__name__}: {e}")
    print(f"\n{len(tests) - failed} passed, {failed} failed")
    plot_results()
    raise SystemExit(failed > 0)
