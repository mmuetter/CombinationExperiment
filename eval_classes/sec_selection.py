import matplotlib

matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os
import sys


class SelectSection:
    def __init__(self, path):
        self.path = path
        self.section_names = os.listdir(path)
        self.clicked_image_name = None  # Attribute to store the clicked image name
        self.should_stop = False  # Flag to control the flow based on user input

    def get_sec_names(self, t, p, sec):
        identifier = f"{t}_{p}"
        sections = [
            name for name in self.section_names if identifier in name and sec in name
        ]
        return sections

    def display_and_pick(self, t, p, sec):
        if self.should_stop:  # Check the flag at the beginning
            return None

        sections = self.get_sec_names(t, p, sec)
        d_values = sorted({name.split("_")[3][1] for name in sections})
        n_values = sorted({name.split("_")[4][1] for name in sections}, key=int)

        fig, axs = plt.subplots(
            len(d_values), len(n_values), squeeze=False, figsize=(18, 12)
        )
        fig.suptitle(t + "_" + p + "_" + sec, fontsize=16)
        image_map = {}  # Map to store axs to image name mapping

        for name in sections:
            parts = name.split("_")
            d_index = d_values.index(parts[3][1])
            n_index = n_values.index(parts[4][1])
            img_path = os.path.join(self.path, name)
            img = mpimg.imread(img_path)
            axs[d_index, n_index].imshow(img)
            axs[d_index, n_index].axis("off")
            axs[d_index, n_index].set_title(f"d={parts[3][1]}, n={parts[4][1]}")
            image_map[axs[d_index, n_index]] = name  # Store reference to image name

        def onclick(event):
            ax = event.inaxes
            if ax and ax in image_map:
                self.clicked_image_name = image_map[ax]
                plt.close(fig)

        fig.canvas.mpl_connect("button_press_event", onclick)

        plt.show()
        return self.clicked_image_name
