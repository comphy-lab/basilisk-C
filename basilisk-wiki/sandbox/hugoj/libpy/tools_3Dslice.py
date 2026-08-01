import numpy as np
import pyvista as pv


def index_slice(grid, axis, index, size):
    """
    Extract a single-index slab of the grid along one axis.

    Inputs
    ----------
    grid : pv.StructuredGrid
        Grid with dims (nx+1, ny+1, nz+1) points, as built by build_3Dgrid.
    axis : {"x", "y", "z"}
    index : int
        Point index along the chosen axis (0..size[axis]).
    size : (nx, ny, nz)
        Cell counts, same convention as surface_extractor/add_mesh_domain_side.

    Returns
    -------
    pv.UnstructuredGrid
        Single-layer slab with cell_data carried over.
    """
    nx, ny, nz = size
    index = int(np.clip(index, 0, {"x": nx, "y": ny, "z": nz}[axis]))
    if axis == "x":
        bounds = (index, index, 0, ny, 0, nz)
    elif axis == "y":
        bounds = (0, nx, index, index, 0, nz)
    elif axis == "z":
        bounds = (0, nx, 0, ny, index, index)
    else:
        raise ValueError(f"axis must be 'x', 'y' or 'z', got {axis!r}")
    return grid.extract_subset(bounds)


class SliceExplorer:
    """
    Adds three axis-select buttons + a slider to a pv.Plotter, letting you
    interactively pick out an exact layer/column slice of `grid` and view
    a scalar field on it.
    """

    def __init__(self, plotter, grid, size, var, start_axis="z", start_index=None):
        """
        plotter : pv.Plotter
        grid : pv.StructuredGrid
        size : (nx, ny, nz)
        var : dict
            {"name": ..., "clim": [...], "cmap": ...} - the scalar field to show.
        """
        self.plotter = plotter
        self.grid = grid
        self.nx, self.ny, self.nz = size
        self.var = var

        self.axis = start_axis
        self.index = self._max_index(self.axis) if start_index is None else start_index

        self.actor = None
        self.buttons = {}
        self.slider = None

        self._build_axis_buttons()
        self._build_slider()
        self._update_slice()

    # -- helpers -----------------------------------------------------
    def _max_index(self, axis):
        return {"x": self.nx, "y": self.ny, "z": self.nz}[axis]

    # -- axis buttons --------------------------------------------------
    def _build_axis_buttons(self):
        positions = {"x": (10, 10), "y": (10, 50), "z": (10, 90)}
        for ax, pos in positions.items():
            self.plotter.add_text(
                ax.upper(),
                position=(pos[0] + 40, pos[1]),
                font_size=10,
                name=f"label_{ax}",
            )
            btn = self.plotter.add_checkbox_button_widget(
                callback=lambda state, axis=ax: self._on_axis_button(axis, state),
                value=(ax == self.axis),
                position=pos,
                size=30,
                color_on="lightblue",
                color_off="grey",
            )
            self.buttons[ax] = btn

    def _on_axis_button(self, axis, state):
        # Treat the three checkboxes as a radio group: clicking the
        # already-active one just snaps it back on; clicking another
        # switches axis and turns the rest off.
        if not state:
            self.buttons[axis].GetRepresentation().SetState(1)
            return
        if axis == self.axis:
            return
        self.axis = axis
        for other_ax, btn in self.buttons.items():
            if other_ax != axis:
                btn.GetRepresentation().SetState(0)
        self.index = min(self.index, self._max_index(axis))
        self._rebuild_slider()
        self._update_slice()

    # -- slider ----------------------------------------------------------
    def _build_slider(self):
        self.slider = self.plotter.add_slider_widget(
            callback=self._on_slider,
            rng=[0, self._max_index(self.axis)],
            value=self.index,
            title=f"{self.axis}-index",
            pointa=(0.25, 0.92),
            pointb=(0.75, 0.92),
            style="modern",
            fmt="%.0f",
        )

    def _rebuild_slider(self):
        # PyVista has no "update range" API for an existing slider,
        # so drop it and add a fresh one sized for the new axis.
        self.plotter.clear_slider_widgets()
        self._build_slider()

    def _on_slider(self, value):
        self.index = int(round(value))
        self._update_slice()

    # -- drawing -----------------------------------------------------
    def _update_slice(self):
        sl = index_slice(self.grid, self.axis, self.index, (self.nx, self.ny, self.nz))
        if self.actor is None:
            self.actor = self.plotter.add_mesh(
                sl,
                scalars=self.var["name"],
                cmap=self.var["cmap"],
                clim=self.var["clim"],
                scalar_bar_args={"title": self.var["name"], "fmt": "%.3f"},
            )
        else:
            # Swap the data in place instead of remove+add, so the
            # scalar bar/camera don't jump on every slider tick.
            self.actor.mapper.SetInputData(sl)
            self.plotter.render()


# Maps each spatial axis to the index-variable name used in the label,
# e.g. "z= -52.4 (l=10)"
_INDEX_NAME = {"x": "i", "y": "j", "z": "l"}


class SliceExplorer2:
    """
    Adds three axis-select buttons, a pair of -/+ index buttons, and a
    text readout to a pv.Plotter, letting you interactively pick out an
    exact layer/column slice of `grid` and view a scalar field on it.
    """

    def __init__(
        self, plotter, grid, size, var, xv, yv, zv, start_axis="z", start_index=None
    ):
        """
        plotter : pv.Plotter
        grid : pv.StructuredGrid
        size : (nx, ny, nz)
        var : dict
            {"name": ..., "clim": [...], "cmap": ...} - the scalar field to show.
        xv, yv : (nx+1,), (ny+1,)
            Vertex coordinates along x and y (from vertex_coordinates).
        zv : (nz+1, ny+1, nx+1)
            Vertex elevations (from interfaces_to_vertices). Since the grid
            is topography-following, elevation varies across a level, so
            the z position shown is the mean over that level's vertices.
        """
        self.plotter = plotter
        self.grid = grid
        self.nx, self.ny, self.nz = size
        self.var = var
        self.xv, self.yv, self.zv = xv, yv, zv

        self.axis = start_axis
        self.index = self._max_index(self.axis) if start_index is None else start_index

        self.actor = None
        self.axis_buttons = {}
        self.pm_buttons = {}

        self._build_axis_buttons()
        self._build_index_buttons()
        self._update_slice()

    # -- helpers -----------------------------------------------------
    def _max_index(self, axis):
        return {"x": self.nx, "y": self.ny, "z": self.nz}[axis]

    def _position_value(self):
        if self.axis == "x":
            return self.xv[self.index]
        elif self.axis == "y":
            return self.yv[self.index]
        else:
            return float(np.mean(self.zv[self.index]))

    # -- axis buttons (radio group) --------------------------------------
    def _build_axis_buttons(self):
        positions = {"x": (10, 10), "y": (10, 50), "z": (10, 90)}
        for ax, pos in positions.items():
            self.plotter.add_text(
                ax.upper(),
                position=(pos[0] + 40, pos[1]),
                font_size=10,
                name=f"label_{ax}",
            )
            btn = self.plotter.add_checkbox_button_widget(
                callback=lambda state, axis=ax: self._on_axis_button(axis, state),
                value=(ax == self.axis),
                position=pos,
                size=30,
                color_on="lightblue",
                color_off="grey",
            )
            self.axis_buttons[ax] = btn

    def _on_axis_button(self, axis, state):
        # Treat the three checkboxes as a radio group: clicking the
        # already-active one just snaps it back on; clicking another
        # switches axis and turns the rest off.
        if not state:
            self.axis_buttons[axis].GetRepresentation().SetState(1)
            return
        if axis == self.axis:
            return
        self.axis = axis
        for other_ax, btn in self.axis_buttons.items():
            if other_ax != axis:
                btn.GetRepresentation().SetState(0)
        self.index = min(self.index, self._max_index(axis))
        self._update_slice()

    # -- index +/- buttons -------------------------------------------
    def _build_index_buttons(self):
        # "-" button
        self.plotter.add_text("-", position=(10, 140), font_size=14, name="label_minus")
        self.pm_buttons["-"] = self.plotter.add_checkbox_button_widget(
            callback=self._on_minus,
            value=False,
            position=(10, 130),
            size=30,
            color_on="lightgreen",
            color_off="grey",
        )
        # "+" button
        self.plotter.add_text("+", position=(60, 140), font_size=14, name="label_plus")
        self.pm_buttons["+"] = self.plotter.add_checkbox_button_widget(
            callback=self._on_plus,
            value=False,
            position=(60, 130),
            size=30,
            color_on="lightgreen",
            color_off="grey",
        )

    def _on_minus(self, state):
        if state:
            self.index = max(0, self.index - 1)
            self._update_slice()
        # momentary push: snap back off regardless of index change
        self.pm_buttons["-"].GetRepresentation().SetState(0)

    def _on_plus(self, state):
        if state:
            self.index = min(self._max_index(self.axis), self.index + 1)
            self._update_slice()
        self.pm_buttons["+"].GetRepresentation().SetState(0)

    # -- drawing -----------------------------------------------------
    def _update_label(self):
        val = self._position_value()
        idx_name = _INDEX_NAME[self.axis]
        text = f"{self.axis}= {val:.1f} ({idx_name}={self.index})"
        self.plotter.add_text(
            text, position=(200, 130), font_size=14, name="slice_position_label"
        )

    def _update_slice(self):
        sl = index_slice(self.grid, self.axis, self.index, (self.nx, self.ny, self.nz))
        if self.actor is None:
            self.actor = self.plotter.add_mesh(
                sl,
                scalars=self.var["name"],
                cmap=self.var["cmap"],
                clim=self.var["clim"],
                scalar_bar_args={"title": self.var["name"], "fmt": "%.3f"},
            )
        else:
            # Swap the data in place instead of remove+add, so the
            # scalar bar/camera don't jump on every click.
            self.actor.mapper.SetInputData(sl)
            self.plotter.render()
        self._update_label()
