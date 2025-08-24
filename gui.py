"""
Main GUI module for the Digital Exoplanet Hunting Observatory.
"""

import tkinter as tk
from tkinter import filedialog, messagebox
import tkinter.ttk as ttk
import logging
import threading
import os

try:
    from PIL import Image, ImageTk
    PIL_AVAILABLE = True
except ImportError:
    PIL_AVAILABLE = False
    print("PIL not available - image display will be limited")

from config import APP_CONFIG, FILE_TYPES, DEFAULT_PARAMS
from photometry import Photometry
from conversion import Conversion
from fits_handler import FitsHandler
from database import ObservatoryDatabase

logger = logging.getLogger(__name__)


class MainGUI:
    """
    Main GUI class for the Digital Exoplanet Hunting Observatory application.
    """
    
    def __init__(self):
        """Initialize the main GUI application."""
        self.root = None
        self.statusbar = None
        self.progressbar = None
        self.progressbar2 = None
        self.progress_window = None
        self.photometry_window = None
        
        # Initialize modules
        self.photometry = Photometry(gui_reference=self)
        self.conversion = Conversion(gui_reference=self)
        self.fits_handler = FitsHandler()
        
        # File paths
        self.current_image_path = None
        self.current_fits_path = None
        self.current_raw_path = None
        self.fits_files = []
        
        logger.info("GUI initialized")
    
    def create_main_window(self):
        """Create and configure the main application window."""
        self.root = tk.Tk()
        self.root.geometry(f"{APP_CONFIG['window_size']}{APP_CONFIG['window_position']}")
        self.root.title(APP_CONFIG['name'])
        self.root.resizable(True, True)
        
        # Set window icon (if available)
        try:
            self.root.iconbitmap('icon.ico')
        except:
            pass  # Icon file not found, continue without it
        
        # Create menu bar
        self._create_menu_bar()
        
        # Create main content area
        self._create_main_content()
        
        # Create status bar
        self._create_status_bar()
        
        logger.info("Main window created")
    
    def _create_menu_bar(self):
        """Create the application menu bar."""
        main_menu = tk.Menu(self.root)
        self.root.config(menu=main_menu)
        
        # File Menu
        file_menu = tk.Menu(main_menu, tearoff=0)
        main_menu.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="Open Image", command=self.open_image, accelerator="Ctrl+O")
        file_menu.add_command(label="Open RAW File", command=self.open_raw_file)
        file_menu.add_command(label="Open FITS File", command=self.open_fits_file)
        file_menu.add_separator()
        file_menu.add_command(label="Show Current Image", command=self.show_current_image, state=tk.DISABLED)
        file_menu.add_command(label="Show FITS Data", command=self.show_fits_data, state=tk.DISABLED)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.exit_application, accelerator="Ctrl+Q")
        
        # Process Menu
        process_menu = tk.Menu(main_menu, tearoff=0)
        main_menu.add_cascade(label="Process", menu=process_menu)
        process_menu.add_command(label="Convert RAW to FITS", command=self.convert_raw_to_fits)
        process_menu.add_command(label="Stack FITS Images", command=self.stack_fits_images)
        
        # Analysis Menu
        analysis_menu = tk.Menu(main_menu, tearoff=0)
        main_menu.add_cascade(label="Analysis", menu=analysis_menu)
        analysis_menu.add_command(label="Aperture Photometry", command=self.start_photometry)
        
        # Database Menu
        database_menu = tk.Menu(main_menu, tearoff=0)
        main_menu.add_cascade(label="Database", menu=database_menu)
        database_menu.add_command(label="Clear All Database Entries", command=self.clear_database)
        
        # Help Menu
        help_menu = tk.Menu(main_menu, tearoff=0)
        main_menu.add_cascade(label="Help", menu=help_menu)
        help_menu.add_command(label="User Guide", command=self.show_help)
        help_menu.add_separator()
        help_menu.add_command(label="About", command=self.show_about)
        
        # Store menu references for enabling/disabling
        self.file_menu = file_menu
    
    def _create_main_content(self):
        """Create the main content area of the application."""
        # Create main frame
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Welcome label
        welcome_label = ttk.Label(
            main_frame, 
            text=f"Welcome to {APP_CONFIG['name']}", 
            font=('Arial', 14, 'bold')
        )
        welcome_label.pack(pady=20)
        
        # Quick action buttons frame
        buttons_frame = ttk.LabelFrame(main_frame, text="Quick Actions")
        buttons_frame.pack(fill=tk.X, pady=10)
        
        # Create button grid
        button_data = [
            ("Open Image", self.open_image, "Open an image file for viewing"),
            ("Open RAW", self.open_raw_file, "Open a RAW astronomical image"),
            ("Open FITS", self.open_fits_file, "Open a FITS astronomical image"),
            ("Convert to FITS", self.convert_raw_to_fits, "Convert RAW image to FITS format"),
            ("Stack Images", self.stack_fits_images, "Combine multiple FITS images"),
            ("Photometry", self.start_photometry, "Perform stellar photometry analysis")
        ]
        
        for i, (text, command, tooltip) in enumerate(button_data):
            btn = ttk.Button(buttons_frame, text=text, command=command, width=20)
            row = i // 3
            col = i % 3
            btn.grid(row=row, column=col, padx=5, pady=5, sticky='ew')
            
            # Add tooltip (simple version)
            self._create_tooltip(btn, tooltip)
        
        # Configure grid weights
        for i in range(3):
            buttons_frame.columnconfigure(i, weight=1)
        
        # Information frame
        info_frame = ttk.LabelFrame(main_frame, text="Current Session")
        info_frame.pack(fill=tk.BOTH, expand=True, pady=10)
        
        self.info_text = tk.Text(info_frame, height=8, state=tk.DISABLED, wrap=tk.WORD)
        scrollbar = ttk.Scrollbar(info_frame, orient=tk.VERTICAL, command=self.info_text.yview)
        self.info_text.configure(yscrollcommand=scrollbar.set)
        
        self.info_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        # Add initial message
        self._add_info_message("Application started. Ready for astronomical image processing.")
    
    def _create_status_bar(self):
        """Create the status bar at the bottom of the window."""
        self.statusbar = ttk.Label(
            self.root, 
            text="Ready", 
            relief=tk.SUNKEN, 
            anchor=tk.W,
            padding=5
        )
        self.statusbar.pack(side=tk.BOTTOM, fill=tk.X)
    
    def _create_tooltip(self, widget, text):
        """Create a simple tooltip for a widget."""
        def show_tooltip(event):
            tooltip = tk.Toplevel()
            tooltip.wm_overrideredirect(True)
            tooltip.wm_geometry(f"+{event.x_root+10}+{event.y_root+10}")
            label = ttk.Label(tooltip, text=text, background="lightyellow", relief=tk.SOLID)
            label.pack()
            widget.tooltip = tooltip
        
        def hide_tooltip(event):
            if hasattr(widget, 'tooltip'):
                widget.tooltip.destroy()
                del widget.tooltip
        
        widget.bind("<Enter>", show_tooltip)
        widget.bind("<Leave>", hide_tooltip)
    
    def _add_info_message(self, message):
        """Add a message to the information text area."""
        self.info_text.config(state=tk.NORMAL)
        self.info_text.insert(tk.END, f"• {message}\\n")
        self.info_text.see(tk.END)
        self.info_text.config(state=tk.DISABLED)
        logger.info(message)
    
    def _update_status(self, message):
        """Update the status bar message."""
        if self.statusbar:
            self.statusbar.config(text=message)
            self.root.update_idletasks()
    
    # File operations
    def open_image(self):
        """Open an image file for viewing."""
        try:
            filepath = filedialog.askopenfilename(
                title="Select Image File",
                filetypes=FILE_TYPES['images']
            )
            
            if filepath:
                self.current_image_path = filepath
                self.file_menu.entryconfig("Show Current Image", state=tk.NORMAL)
                self._add_info_message(f"Image loaded: {os.path.basename(filepath)}")
                self._update_status(f"Image loaded: {os.path.basename(filepath)}")
                
        except Exception as e:
            logger.error(f"Error opening image: {str(e)}")
            messagebox.showerror('Error', f'Failed to open image: {e}')
    
    def open_raw_file(self):
        """Open a RAW file."""
        try:
            filepath = filedialog.askopenfilename(
                title="Select RAW File",
                filetypes=FILE_TYPES['raw']
            )
            
            if filepath:
                self.current_raw_path = filepath
                self.current_image_path = filepath  # Enable RAW files to be shown as images
                self.file_menu.entryconfig("Show Current Image", state=tk.NORMAL)
                self._add_info_message(f"RAW file loaded: {os.path.basename(filepath)}")
                self._update_status(f"RAW file ready for conversion")
                
        except Exception as e:
            logger.error(f"Error opening RAW file: {str(e)}")
            messagebox.showerror('Error', f'Failed to open RAW file: {e}')
    
    def open_fits_file(self):
        """Open a FITS file."""
        try:
            filepath = filedialog.askopenfilename(
                title="Select FITS File",
                filetypes=FILE_TYPES['fits']
            )
            
            if filepath:
                self.current_fits_path = filepath
                self.file_menu.entryconfig("Show FITS Data", state=tk.NORMAL)
                self._add_info_message(f"FITS file loaded: {os.path.basename(filepath)}")
                self._update_status(f"FITS file ready for analysis")
                
        except Exception as e:
            logger.error(f"Error opening FITS file: {str(e)}")
            messagebox.showerror('Error', f'Failed to open FITS file: {e}')
    
    def show_current_image(self):
        """Display the currently loaded image."""
        if not self.current_image_path:
            messagebox.showwarning('No Image', 'No image file is currently loaded.')
            return
        
        try:
            self._show_image_in_window(self.current_image_path)
        except Exception as e:
            logger.error(f"Error displaying image: {str(e)}")
            messagebox.showerror('Display Error', f'Failed to display image: {e}')
    
    def show_fits_data(self):
        """Display the currently loaded FITS file."""
        if not self.current_fits_path:
            messagebox.showwarning('No FITS File', 'No FITS file is currently loaded.')
            return
        
        try:
            self.fits_handler.load_fits(self.current_fits_path)
            self._add_info_message(f"Displayed FITS data: {os.path.basename(self.current_fits_path)}")
        except Exception as e:
            logger.error(f"Error displaying FITS: {str(e)}")
            messagebox.showerror('Display Error', f'Failed to display FITS: {e}')
    
    def _show_image_in_window(self, filepath):
        """Show image in a separate window with scrollbars."""
        if not PIL_AVAILABLE:
            messagebox.showerror('PIL Not Available', 'PIL/Pillow is required for image display. Please install it with: pip install Pillow')
            return
        
        try:
            img_window = tk.Toplevel(self.root)
            img_window.title(f"Image: {os.path.basename(filepath)}")
            
            # Create scrollable canvas
            canvas = tk.Canvas(img_window, bg='black')
            v_scrollbar = ttk.Scrollbar(img_window, orient=tk.VERTICAL, command=canvas.yview)
            h_scrollbar = ttk.Scrollbar(img_window, orient=tk.HORIZONTAL, command=canvas.xview)
            canvas.configure(yscrollcommand=v_scrollbar.set, xscrollcommand=h_scrollbar.set)
            
            # Load and display image
            pil_image = Image.open(filepath)
            # Resize if too large
            max_size = (1200, 900)
            pil_image.thumbnail(max_size, Image.Resampling.LANCZOS)
            
            photo = ImageTk.PhotoImage(pil_image)
            canvas.create_image(0, 0, anchor=tk.NW, image=photo)
            canvas.configure(scrollregion=canvas.bbox("all"))
            
            # Pack elements
            canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
            v_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
            h_scrollbar.pack(side=tk.BOTTOM, fill=tk.X)
            
            # Keep a reference to prevent garbage collection
            canvas.image = photo
            
            img_window.geometry("800x600")
            
        except Exception as e:
            logger.error(f"Error displaying image: {str(e)}")
            messagebox.showerror('Display Error', f'Failed to display image: {e}')
    
    # Processing operations
    def convert_raw_to_fits(self):
        """Convert RAW file to FITS format."""
        if not self.current_raw_path:
            # Ask user to select RAW file
            self.open_raw_file()
            if not self.current_raw_path:
                return
        
        # Ask for save location
        save_path = filedialog.asksaveasfilename(
            title="Save FITS File",
            defaultextension=".fits",
            filetypes=FILE_TYPES['fits']
        )
        
        if save_path:
            self._show_progress_window("RAW to FITS Conversion")
            
            # Run conversion in a separate thread to prevent GUI freezing
            def conversion_thread():
                success = self.conversion.raw_to_fits(self.current_raw_path, save_path)
                if success:
                    self._add_info_message(f"RAW converted to FITS: {os.path.basename(save_path)}")
            
            threading.Thread(target=conversion_thread, daemon=True).start()
    
    def stack_fits_images(self):
        """Stack multiple FITS images."""
        filepaths = filedialog.askopenfilenames(
            title="Select FITS Files to Stack",
            filetypes=FILE_TYPES['fits']
        )
        
        if filepaths:
            self.fits_files = filepaths
            
            # Run stacking in a separate thread
            def stacking_thread():
                result = self.fits_handler.mean_fits(filepaths)
                if result is not None:
                    self._add_info_message(f"Stacked {len(filepaths)} FITS images")
                    
                    # Ask if user wants to save the result
                    if messagebox.askyesno("Save Stacked Image", "Would you like to save the stacked image?"):
                        save_path = filedialog.asksaveasfilename(
                            title="Save Stacked FITS",
                            defaultextension=".fits",
                            filetypes=FILE_TYPES['fits']
                        )
                        if save_path:
                            self.fits_handler.save_fits(result, save_path)
            
            threading.Thread(target=stacking_thread, daemon=True).start()
    
    def start_photometry(self):
        """Start aperture photometry analysis."""
        if not self.current_fits_path:
            # Ask user to select FITS file
            self.open_fits_file()
            if not self.current_fits_path:
                return
        
        self._show_photometry_dialog()
    
    def _show_photometry_dialog(self):
        """Show photometry parameters dialog."""
        self.photometry_window = tk.Toplevel(self.root)
        self.photometry_window.title("Aperture Photometry Parameters")
        self.photometry_window.geometry("350x200")
        self.photometry_window.resizable(False, False)
        
        # Make window modal
        self.photometry_window.transient(self.root)
        self.photometry_window.grab_set()
        
        # FWHM parameter
        ttk.Label(self.photometry_window, text="FWHM (Full Width Half Maximum):").pack(pady=5)
        self.fwhm_var = tk.StringVar(value=str(DEFAULT_PARAMS['fwhm']))
        fwhm_entry = ttk.Entry(self.photometry_window, textvariable=self.fwhm_var, width=25)
        fwhm_entry.pack(pady=5)
        
        # Threshold parameter
        ttk.Label(self.photometry_window, text="Detection Threshold (sigma):").pack(pady=5)
        self.threshold_var = tk.StringVar(value=str(DEFAULT_PARAMS['threshold']))
        threshold_entry = ttk.Entry(self.photometry_window, textvariable=self.threshold_var, width=25)
        threshold_entry.pack(pady=5)
        
        # Buttons frame
        buttons_frame = ttk.Frame(self.photometry_window)
        buttons_frame.pack(pady=20)
        
        ttk.Button(buttons_frame, text="Start Analysis", command=self._perform_photometry).pack(side=tk.LEFT, padx=5)
        ttk.Button(buttons_frame, text="Cancel", command=self.photometry_window.destroy).pack(side=tk.LEFT, padx=5)
        
        # Focus on first entry
        fwhm_entry.focus()
    
    def _perform_photometry(self):
        """Perform the photometry analysis with user parameters."""
        try:
            fwhm = float(self.fwhm_var.get())
            threshold = float(self.threshold_var.get())
            
            if fwhm <= 0 or threshold <= 0:
                messagebox.showerror('Invalid Parameters', 'FWHM and threshold must be positive values.')
                return
            
            self.photometry_window.destroy()
            self._show_photometry_progress_window()
            
            # Run photometry in separate thread
            def photometry_thread():
                sources, phot_table = self.photometry.photometry(self.current_fits_path, fwhm, threshold)
                if sources is not None:
                    self._add_info_message(f"Photometry completed: {len(sources)} sources found")
            
            threading.Thread(target=photometry_thread, daemon=True).start()
        except ValueError:
            messagebox.showerror('Invalid Input', 'Please enter valid numeric values for FWHM and threshold.')
        except Exception as e:
            logger.error(f"Error in photometry: {str(e)}")
            messagebox.showerror('Photometry Error', f'Failed to perform photometry: {e}')
            
    def clear_database(self):
        """Clear all entries from the database after user confirmation."""
        # Show a very explicit warning to the user
        confirm = messagebox.askyesno(
            "Confirm Database Deletion",
            "WARNING: This will permanently delete ALL data from the database, "
            "including all sessions, sources, and catalog matches.\n\n"
            "This action cannot be undone.\n\n"
            "Are you absolutely sure you want to continue?",
            icon='warning',
            default=messagebox.NO
        )

        if confirm:
            logger.warning("User confirmed database deletion.")
            try:
                db = ObservatoryDatabase()
                if not db.initialize_connection():
                    messagebox.showerror("Database Error", "Could not connect to the database.")
                    self._add_info_message("Failed to connect to database for clearing.")
                    return

                success = db.clear_all_data()
                db.close()

                if success:
                    messagebox.showinfo("Database Cleared", "All entries have been successfully removed from the database.")
                    self._add_info_message("Database has been cleared by user.")
                else:
                    messagebox.showerror("Database Error", "An error occurred while clearing the database. Please check the logs for details.")
                    self._add_info_message("Error occurred while clearing the database.")

            except Exception as e:
                logger.error(f"Failed to clear database from GUI: {e}", exc_info=True)
                messagebox.showerror("Database Error", f"A critical error occurred while clearing the database:\n\n{e}")
                self._add_info_message(f"Critical error during database clearing: {e}")
    
    def _show_progress_window(self, title):
        """Show a progress window for long operations."""
        self.progress_window = tk.Toplevel(self.root)
        self.progress_window.title(title)
        self.progress_window.geometry("400x100")
        self.progress_window.resizable(False, False)
        
        # Make window modal
        self.progress_window.transient(self.root)
        self.progress_window.grab_set()
        
        # Progress bar
        ttk.Label(self.progress_window, text=f"{title} in progress...").pack(pady=10)
        self.progressbar = ttk.Progressbar(
            self.progress_window, 
            mode='determinate', 
            length=350
        )
        self.progressbar.pack(pady=10)
        
        # Center the window
        self.progress_window.update_idletasks()
        x = (self.progress_window.winfo_screenwidth() // 2) - (self.progress_window.winfo_width() // 2)
        y = (self.progress_window.winfo_screenheight() // 2) - (self.progress_window.winfo_height() // 2)
        self.progress_window.geometry(f"+{x}+{y}")
    
    def _show_photometry_progress_window(self):
        """Show progress window for photometry analysis."""
        self.progress_window = tk.Toplevel(self.root)
        self.progress_window.title("Photometry Analysis")
        self.progress_window.geometry("400x100")
        self.progress_window.resizable(False, False)
        
        # Make window modal
        self.progress_window.transient(self.root)
        self.progress_window.grab_set()
        
        # Progress bar
        ttk.Label(self.progress_window, text="Performing aperture photometry...").pack(pady=10)
        self.progressbar2 = ttk.Progressbar(
            self.progress_window, 
            mode='determinate', 
            length=350
        )
        self.progressbar2.pack(pady=10)
        
        # Store reference for cleanup
        self.phot_barwindow = self.progress_window
        
        # Center the window
        self.progress_window.update_idletasks()
        x = (self.progress_window.winfo_screenwidth() // 2) - (self.progress_window.winfo_width() // 2)
        y = (self.progress_window.winfo_screenheight() // 2) - (self.progress_window.winfo_height() // 2)
        self.progress_window.geometry(f"+{x}+{y}")
    
    # Help and About
    def show_help(self):
        """Show user guide/help information."""
        help_window = tk.Toplevel(self.root)
        help_window.title("User Guide")
        help_window.geometry("600x400")
        
        help_text = tk.Text(help_window, wrap=tk.WORD, padx=10, pady=10)
        scrollbar = ttk.Scrollbar(help_window, orient=tk.VERTICAL, command=help_text.yview)
        help_text.configure(yscrollcommand=scrollbar.set)
        
        help_content = """
Digital Exoplanet Hunting Observatory - User Guide

1. OPENING FILES
   • Use File menu to open images, RAW files, or FITS files
   • Supported formats: PNG, TIFF, JPEG, GIF, CR2, NEF, ARW, DNG, FITS

2. RAW TO FITS CONVERSION
   • Open a RAW file (File > Open RAW File)
   • Select Process > Convert RAW to FITS
   • Choose save location for the FITS file
   • Metadata from RAW file is preserved

3. IMAGE STACKING
   • Select Process > Stack FITS Images
   • Choose multiple FITS files to combine
   • The mean-combined result reduces noise

4. APERTURE PHOTOMETRY
   • Open a FITS file for analysis
   • Select Analysis > Aperture Photometry
   • Set FWHM (typical: 3-8 pixels)
   • Set threshold (typical: 2-5 sigma)
   • Results are saved as CSV files

5. TIPS
   • Use FITS format for scientific analysis
   • Stack multiple images to improve signal-to-noise ratio
   • Adjust photometry parameters based on your image quality
   • Check the information panel for processing updates

For technical support, refer to the documentation or contact the development team.
        """
        
        help_text.insert(tk.END, help_content)
        help_text.config(state=tk.DISABLED)
        
        help_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
    
    def show_about(self):
        """Show about dialog."""
        about_window = tk.Toplevel(self.root)
        about_window.title("About")
        about_window.geometry("400x250")
        about_window.resizable(False, False)
        
        # Center the window
        about_window.update_idletasks()
        x = (about_window.winfo_screenwidth() // 2) - (about_window.winfo_width() // 2)
        y = (about_window.winfo_screenheight() // 2) - (about_window.winfo_height() // 2)
        about_window.geometry(f"+{x}+{y}")
        
        # Content
        ttk.Label(about_window, text=APP_CONFIG['about_text']['name'], 
                 font=('Arial', 14, 'bold')).pack(pady=10)
        ttk.Label(about_window, text=f"Version {APP_CONFIG['about_text']['version']}").pack()
        
        ttk.Label(about_window, text="Developed by:", font=('Arial', 10, 'bold')).pack(pady=(20, 5))
        for dev in APP_CONFIG['about_text']['developers']:
            ttk.Label(about_window, text=f"• {dev}").pack()
        
        ttk.Label(about_window, text=APP_CONFIG['about_text']['institution'], 
                 font=('Arial', 9)).pack(pady=20)
        
        ttk.Button(about_window, text="Close", command=about_window.destroy).pack(pady=10)
    
    def exit_application(self):
        """Exit the application."""
        if messagebox.askyesno("Exit", "Are you sure you want to exit?"):
            logger.info("Application closing")
            self.root.destroy()
    
    def run(self):
        """Start the GUI application."""
        self.create_main_window()
        
        # Bind keyboard shortcuts
        self.root.bind('<Control-o>', lambda e: self.open_image())
        self.root.bind('<Control-q>', lambda e: self.exit_application())
        
        # Handle window closing
        self.root.protocol("WM_DELETE_WINDOW", self.exit_application)
        
        logger.info("Starting GUI main loop")
        self.root.mainloop()