from __future__ import annotations
import customtkinter as ctk
from app.ui.tab_view import PeptideTabView


class App(ctk.CTk):
    """Main application window for the Peptide Sequence Tool."""

    def __init__(self) -> None:
        super().__init__()
        self.title("Peptide Sequence Tool")
        ctk.set_appearance_mode("dark")

        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()
        self.geometry(f"{screen_width}x{screen_height}")

        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)

        self.output_text = ctk.CTkTextbox(self, height=360, width=480)
        self.output_text.grid(row=1, column=0, padx=10, pady=10, sticky="nsew")

        self.tabview = PeptideTabView(self, self.output_text)
        self.tabview.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")