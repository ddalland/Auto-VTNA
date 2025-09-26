import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from auto_vtna.VTNA_functions import same_excess_v2

class Same_Excess:
    def __init__(self, data, selected_experiments, selected_species, selected_product_species=None):
        """
        Parameters
        ----------
        data : dict[str, pd.DataFrame]
            Each DataFrame must have time as its first column, then species columns.
        selected_experiments : list[str]
            Keys in `data` to include.
        selected_species : str
            Reactant to use for Same-Excess alignment (must decrease overall).
        """
        self.data = data
        self.selected_experiments = list(selected_experiments)
        self.selected_species = selected_species
        self.selected_product_species=selected_product_species

        # Run the generalized same_excess (works with 2+ experiments).
        results = same_excess_v2(self.data, self.selected_experiments, self.selected_species,selected_product_species=self.selected_product_species)
        if isinstance(results, str):
            self.result = None
            self.error = results
            return
        else:
            self.result, self.low_initial_conc_experiment = results 
            self.error = None

    def plot_same_excess(self, y_unit='M', size_scaler=1, title=None, grid=False,
                     xaxis_size_scaler=None, DP_scaler=0.35, xtick_rotation=0,
                     show_legend=True, legend_outside=False, plt_show=True,
                     linewidth=0.7, xaxis_notation='Automatic',
                     title_fontsize=10,lim_xy_0=False):

        plt.rcParams['font.family'] = 'DejaVu Sans'

        if isinstance(self.result, str) or self.result is None:
            print(self.error if self.error else "Same_Excess: no result to plot.")
            return

        fig_size = (5 * size_scaler, 4 * size_scaler)
        fig, ax = plt.subplots(figsize=fig_size)

        if DP_scaler is None or isinstance(DP_scaler, str):
            DP_scaler = 0.35

        example_exp = self.selected_experiments[0]
        time_label = self.result[example_exp].columns[0]

        # Standard experiment: highest initial reactant
        init_vals = {exp: float(self.result[exp][self.selected_species].iloc[0])
                    for exp in self.selected_experiments}
        standard_exp = max(init_vals, key=init_vals.get)

        # choose what to plot on y
        y_col = self.selected_product_species or self.selected_species

        plots, base_zero_time, init_shift_val = {}, {}, {}

        for exp in self.selected_experiments:
            df = self.result[exp]
            x = np.asarray(df[time_label], dtype=float)
            y = np.asarray(df[y_col], dtype=float)

            x0 = float(x[0]) if x.size else 0.0
            base_zero_time[exp] = x - x0
            init_shift_val[exp] = x0

            label = f"{exp} (standard)" if exp == standard_exp else exp
            (line,) = ax.plot(x, y, marker='o', markersize=6 * DP_scaler,
                            linewidth=linewidth, label=label)
            plots[exp] = line

        # --- Sliders ---
        non_standard = [e for e in self.selected_experiments if e != standard_exp]
        n_sliders = len(non_standard)
        sliders = {}

        slider_height = 0.035
        slider_pad = 0.010
        xlabel_gap = 0.08 * size_scaler
        bottom_base = 0.08
        bottom_needed = bottom_base + xlabel_gap + n_sliders * (slider_height + slider_pad)
        bottom_needed = min(0.90, bottom_needed)
        plt.subplots_adjust(bottom=bottom_needed)

        left, width = 0.18, 0.75
        first_slider_bottom = 0.02
        current_bottom = first_slider_bottom

        for exp in non_standard:
            ax_slider = plt.axes([left, current_bottom, width, slider_height],
                                facecolor='lightgoldenrodyellow')

            df = self.result[exp]
            x = np.asarray(df[time_label], dtype=float)
            span = float(x[-1] - x[0]) if len(x) > 1 else max(1.0, abs(init_shift_val[exp]) or 1.0)

            vinit = init_shift_val[exp]
            vmin = vinit - 0.10 * abs(span)
            vmax = vinit + 0.10 * abs(span)
            vstep = abs(span) / 1000.0 if abs(span) > 0 else None

            sliders[exp] = Slider(ax=ax_slider, label=f"Shift: {exp}",
                                valmin=vmin, valmax=vmax, valinit=vinit, valstep=vstep)
            current_bottom += (slider_height + slider_pad)

        def update_all(_):
            for exp in non_standard:
                shift = sliders[exp].val
                plots[exp].set_xdata(base_zero_time[exp] + shift)
            fig.canvas.draw_idle()

        for sl in sliders.values():
            sl.on_changed(update_all)

        # --- Cosmetics ---
        if show_legend:
            if legend_outside:
                # First draw the figure so the legend can report its size accurately
                fig.canvas.draw()

                # Create (temporary) outside legend
                leg = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=True)

                # Measure legend width in display units (pixels)
                renderer = fig.canvas.get_renderer()
                leg_bbox = leg.get_window_extent(renderer=renderer)
                leg_w_px = leg_bbox.width

                # Figure width in pixels
                fig_w_in, fig_h_in = fig.get_size_inches()
                dpi = fig.dpi
                fig_w_px = fig_w_in * dpi

                # Fraction of fig width that the legend needs
                frac = leg_w_px / fig_w_px
                pad_frac = 0.03  # a little breathing room

                # Option A: expand figure width so legend is fully visible
                new_fig_w_in = fig_w_in * (1.0 + frac + pad_frac)
                fig.set_size_inches(new_fig_w_in, fig_h_in, forward=True)

                # Option B (also recommended): tighten the axes so the right side leaves space for the legend
                # Recompute since fig width changed
                fig.canvas.draw()
                # Use the *original* legend width fraction against the *new* fig width to set right margin
                fig.subplots_adjust(right=max(0.5, 1.0 - (frac + pad_frac)))

                # Keep your bottom spacing for sliders and a bit of top margin for the title
                fig.subplots_adjust(bottom=max(0.18, fig.subplotpars.bottom))
                fig.subplots_adjust(top=min(0.85, fig.subplotpars.top))

                # Re-attach legend outside (anchors to the new figure size)
                leg = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=True)
            else:
                ax.legend()


        ax.set_xlabel(time_label, labelpad=7)
        ax.set_ylabel(f'[{y_col}] / {y_unit}')
        ax.set_title(title if title else
                    f"Same excess plot for {y_col}\nStandard: {standard_exp}",
                    fontsize=title_fontsize)
        if lim_xy_0:
            ax.set_ylim(0,None)
            ax.set_xlim(0,None)
        if grid:
            ax.grid()

        if xaxis_notation.lower() in ['scientific', 'sci']:
            ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        elif xaxis_notation.lower() in ['normal']:
            ax.ticklabel_format(style='plain', axis='x')

        ax.tick_params(axis='x', labelrotation=xtick_rotation)
        if not lim_xy_0:
            ax.axhline(0, color='black', linestyle='-', linewidth=0.7)
            ax.axvline(0, color='black', linestyle='-', linewidth=0.7)

        if plt_show:
            plt.show()
        return fig
