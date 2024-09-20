import argparse
import os
import pymol
import re
from pymol import cmd
import subprocess
import shutil
import gemmi
import MDAnalysis as mda
from MDAnalysis.analysis import align
from PIL import Image, ImageDraw, ImageFont
import logging

logger = logging.getLogger(__file__)
logger.setLevel(level=logging.DEBUG)

class ProteinMovieMaker:
    def __init__(self, input_directory, output_movie_file="protein_movie.mp4", output_dcd_file="protein_trajectory.dcd", frame_duration_seconds=1, low_res=False, keep_data=False):
        self.input_directory = input_directory
        self.output_movie_file = output_movie_file
        self.output_dcd_file = output_dcd_file
        self.frame_duration_seconds = frame_duration_seconds
        self.cif_files = sorted(
            [os.path.join(self.input_directory, f) for f in os.listdir(self.input_directory) if f.endswith('.cif')]
        )
        self.pdb_files = []
        self.object_name = "protein_trajectory"
        self.png_dir = os.path.join(input_directory, "temp_png")
        self.pdb_dir = os.path.join(input_directory, "temp_pdb")
        self.low_res = low_res
        self.keep_data = keep_data

    def cif_to_pdb(self):
        os.makedirs(self.pdb_dir, exist_ok=True)

        for cif_file in self.cif_files:
            pdb_file = os.path.join(self.pdb_dir, os.path.basename(cif_file).replace('.cif', '.pdb'))

            # cif to pdb
            doc = gemmi.cif.read_file(cif_file)
            block = doc.sole_block()
            structure = gemmi.make_structure_from_block(block)

            # multiply coordinates by 10 for "evoformer" files
            # TODO this is almost random and should be fixed in openfold and/or dr output...
            if "evoformer" in cif_file:
                for model in structure:
                    for chain in model:
                        for residue in chain:
                            for atom in residue:
                                atom.pos.x *= 10 
                                atom.pos.y *= 10
                                atom.pos.z *= 10

            # structure.remove_hydrogens()  # TODO optional, remove hydrogen if not needed; should we?
            structure.write_pdb(pdb_file)

            logger.debug(f"Successfully converted {cif_file} to {pdb_file}, with scaling applied.")
            self.pdb_files.append(pdb_file)

        if not self.pdb_files:
            raise FileNotFoundError("No pdb files were created from cif files.")



    def center_and_align_pdbs(self):
        centered_pdb_dir = os.path.join(self.pdb_dir, "centered")
        os.makedirs(centered_pdb_dir, exist_ok=True)

        # Use the first structure as reference
        ref_universe = mda.Universe(self.pdb_files[0])
        ref_atoms = ref_universe.select_atoms('name CA')  # I am using alpha carbons for alignment...

        aligned_pdb_files = []

        for pdb_file in self.pdb_files:
            mobile_universe = mda.Universe(pdb_file)
            mobile_atoms = mobile_universe.select_atoms('name CA')

            aligner = align.AlignTraj(
                mobile_universe, ref_universe, select='name CA', in_memory=True
            )
            aligner.run()

            # new pdb file, aligned
            aligned_pdb_file = os.path.join(centered_pdb_dir, os.path.basename(pdb_file))
            with mda.Writer(aligned_pdb_file, multiframe=False) as pdb_writer:
                pdb_writer.write(mobile_universe.atoms)

            aligned_pdb_files.append(aligned_pdb_file)

        self.pdb_files = aligned_pdb_files
        logger.debug(f"Centered and aligned pdb files saved to {centered_pdb_dir}.")

    def load_pdbs(self):
        cmd.reinitialize() 

        # load first pdb file
        cmd.load(self.pdb_files[0], self.object_name)

        # load subsequent pdb files as states
        for i, pdb_file in enumerate(self.pdb_files[1:], start=2):
            cmd.load(pdb_file, self.object_name, state=i)

        logger.debug(f"Loaded {len(self.pdb_files)} pdb files into pymol as states in object '{self.object_name}'.")

    def export_movie(self):
        os.makedirs(self.png_dir, exist_ok=True)

        total_frames = len(self.pdb_files)

        # vis. settings
        
        cmd.hide("everything", self.object_name)
        cmd.dss(self.object_name)  # TODO assign secondary structures; slow, apparently does not work
        # set cartoon to handle loops as fallback if secondary structure fails
        cmd.set("cartoon_loop_quality", 1)  # faster rendering
        cmd.set("cartoon_trace_atoms", 1)  # trace through missing backbone atoms
        cmd.show("cartoon", self.object_name)
        
        cmd.spectrum("count", "rainbow", selection=self.object_name)
        cmd.bg_color("white")
        cmd.set("ray_trace_mode", 0)
        cmd.set("antialias", 1)  #TODO set back to 2 (slower)
        cmd.set("cartoon_transparency", 0.0)
        cmd.set("specular", 0.2)
        cmd.set("ambient", 0.5)
        cmd.set("ray_opaque_background", 1)

        # initial orientation and zoom level
        cmd.orient()
        initial_view = cmd.get_view()  # trying for consistency...

        if self.low_res:
            image_width=800
            image_height=600
        else:        
            image_width = 1920
            image_height = 1080

        x = 0

        for i in range(len(self.pdb_files)):
            state = i + 1
            cmd.frame(state)
            cmd.refresh()

            # restore initial view
            cmd.set_view(initial_view)

            frame_filename = os.path.join(self.png_dir, f"frame{state:04d}.png")
            cmd.png(frame_filename, width=image_width, height=image_height, ray=1)
            logger.info(f"Rendered frame {state}/{total_frames}")


            pdb_file = self.pdb_files[i]
            pdb_basename = os.path.basename(pdb_file)

            # extract "frame" number immediately before _evoformer.pdb or .pdb
            # TODO clearly not the best solution...
            match = re.search(r'_(\d+)(?:_evoformer)?\.pdb$', pdb_basename)
            number = int(match.group(1))

            if "evoformer" in pdb_basename:
                y = number % 48  # 48 blocks in the evoformer stack; let's hope it does not change
                text = f"Recycling iteration {x}, block {y}"
            else: #TODO remove branch to have same output as alphafold
                text = f"Recycling iteration {x}, end"
                x += 1  # next recycle iter

            self.label_image(frame_filename, text)

        logger.debug(f"png frames saved to {self.png_dir}")

        self.create_ffmpeg_list()
        self.pngs_to_mpeg()

    def label_image(self, image_path, text):
        image = Image.open(image_path)
        draw = ImageDraw.Draw(image)

        font_size = 20 if self.low_res else 32

        # load the font from included font file
        font_file = 'LiberationSans-Regular.ttf'
        font_path = os.path.join(os.path.dirname(__file__), font_file)
        if not os.path.exists(font_path):
            print(f"Error: Font file '{font_path}' not found.")
            print(f"Please ensure the font file {font_file} is in the same directory as the script.")
            return

        font = ImageFont.truetype(font_path, font_size)
        max_width = image.width - 40  # 20 pixels padding on each side
        max_height = image.height - 40

        # text in lower-left corner, with padding
        text_position = (20, image.height - font_size - 20)

        draw.text(text_position, text, font=font, fill="black")
        image.save(image_path)

    def create_ffmpeg_list(self):
        list_filename = os.path.join(self.png_dir, 'images.txt')
        with open(list_filename, 'w') as f:
            for i in range(len(self.pdb_files)):
                frame_filename = f"frame{i+1:04d}.png"
                f.write(f"file '{frame_filename}'\n")
                duration_seconds = self.frame_duration_seconds
                f.write(f"duration {duration_seconds}\n")

            # add last image again without duration, because ffmpeg, FFS!!!
            f.write(f"file '{frame_filename}'\n")

        # print(f"Created ffmpeg list file: {list_filename}")

    def pngs_to_mpeg(self):
        list_file_path = os.path.join(self.png_dir, 'images.txt')
        if not os.path.exists(list_file_path):
            print("Error: images.txt file not found in temp_png directory.")
            return

        ffmpeg_cmd = [
            'ffmpeg',
            '-y',  # overwrite output file if exists
            '-f', 'concat',
            '-safe', '0',
            '-i', 'images.txt',
            '-vsync', 'vfr',
            '-pix_fmt', 'yuv420p',
            os.path.abspath(self.output_movie_file)
        ]

        try:
            subprocess.run(ffmpeg_cmd, check=True, cwd=self.png_dir)
            print(f"movie saved as {self.output_movie_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error during movie creation with ffmpeg: {e}")

    def export_traj(self):
        u = mda.Universe(self.pdb_files[0])  # use first pdb as topology
        with mda.Writer(self.output_dcd_file, n_atoms=u.atoms.n_atoms) as dcd_writer:
            for pdb_file in self.pdb_files:
                u = mda.Universe(pdb_file)
                dcd_writer.write(u.atoms)
        print(f"trajectory saved as {self.output_dcd_file}.")

    def clean_up(self):
        if os.path.exists(self.pdb_dir):
            shutil.rmtree(self.pdb_dir)
            print(f"Removed temporary pdb directory: {self.pdb_dir}")
        if os.path.exists(self.png_dir):
            shutil.rmtree(self.png_dir)
            print(f"Removed temporary png directory: {self.png_dir}")

    def run(self):
        # TODO maybe it'd be better having export_movie and export_traj as separate functions:
        # TODO embed first 3 functions just after init...

        # pymol in cmd mode (no gui, quiet)
        pymol.finish_launching(['pymol', '-cq'])
        
        self.cif_to_pdb()
        self.center_and_align_pdbs()
        self.load_pdbs()
        self.export_movie()
        self.export_traj()
        if not self.keep_data:
            self.clean_up()

        # quit pymol
        cmd.quit()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate protein movie and trajectory from cif files")
    parser.add_argument("cif_directory", help="Directory containing cif files")
    parser.add_argument("--output_movie_file", default="protein_movie.mp4", help="Path to the output movie file (default: protein_movie.mp4)")
    parser.add_argument("--output_dcd_file", default="protein_trajectory.dcd", help="Path to the output DCD file (default: protein_trajectory.dcd)")
    parser.add_argument("--frame_duration_seconds", type=int, default=1, help="Frame duration in seconds (default: 1)")
    parser.add_argument("--low_res", action="store_true", help="Generate low resolution movie (default: False)")
    parser.add_argument("--keep_data", action="store_true", help="Keep intermediate data (pdbs and pngs) for debugging (default: False)")

    args = parser.parse_args()

    mmaker = ProteinMovieMaker(
        input_directory=args.cif_directory,
        output_movie_file=args.output_movie_file,
        output_dcd_file=args.output_dcd_file,
        frame_duration_seconds=args.frame_duration_seconds,
        low_res=args.low_res,
        keep_data=args.keep_data
    )
    mmaker.run()

