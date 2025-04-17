import os
import re
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import pandas as pd
from loguru import logger
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas
from reportlab.platypus import (
    BaseDocTemplate,
    Flowable,
    Frame,
    PageTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
)

from rbceq2.core_logic.constants import DB_VERSION, VERSION

# --- Helper Functions ---


def _normalize_sample_id(sample_id: str) -> str:
    """Removes common suffixes from sample IDs for consistent matching.

    Args:
        sample_id: The original sample identifier string.

    Returns:
        The normalized sample identifier string.
    """
    # Ensure input is string before regex
    sample_id_str = str(sample_id)
    # Remove the common long suffix first, including optional .vcf
    normalized = re.sub(
        r"_GRCh38_1_22_v4\.2\.1_benchmark_filtered(\.vcf)?$", "", sample_id_str
    )
    # Example: remove just .vcf if the long suffix wasn't present
    normalized = re.sub(r"\.vcf$", "", normalized)
    # Add any other normalization rules here if needed
    return normalized


def _format_cell_content(text: Optional[Any], separator: Optional[str] = None) -> str:
    """Formats cell text for PDF tables, handling N/A, line breaks, and HTML escaping.

    Args:
        text: The content of the cell (can be string, number, NaN, None).
        separator: The character(s) used to split the text for line breaks (e.g., ',').
                   If None, no splitting occurs, only escaping and N/A handling.

    Returns:
        A string formatted for use in a ReportLab Paragraph, with <br/> for line breaks
        and escaped HTML characters, or "N/A" if the input is missing/empty.
    """
    if pd.isna(text) or text == "" or text is None:
        return "N/A"

    text_str = str(text)
    # Escape HTML special characters first to avoid breaking tags later
    # Using standard library for robust escaping
    import html

    escaped_text = html.escape(text_str)

    if separator:
        # Split the original unescaped text by the separator
        items = text_str.split(separator)
        # Escape *each item* after splitting and strip whitespace
        escaped_items = [html.escape(item.strip()) for item in items]
        # Join the individually escaped items with <br/>
        return "<br/>".join(escaped_items)
    else:
        # Return the fully escaped text if no separator provided
        return escaped_text


# --- Data Preparation ---


def _prepare_dataframes(
    df_genotypes: Optional[pd.DataFrame],
    df_pheno_alpha: Optional[pd.DataFrame],
    df_pheno_num: Optional[pd.DataFrame],
) -> Tuple[Dict[str, Optional[pd.DataFrame]], Set[str], Dict[str, str]]:
    """Normalizes sample IDs in DataFrames and identifies unique samples.

    Adds a 'SampleID_Normalized' column to each input DataFrame.

    Args:
        df_genotypes: DataFrame with genotype data. Index should be SampleID.
        df_pheno_alpha: DataFrame with alphanumeric phenotype data. Index should be SampleID.
        df_pheno_num: DataFrame with numeric phenotype data. Index should be SampleID.

    Returns:
        A tuple containing:
        - A dictionary mapping data type ('genotype', 'alpha', 'numeric') to the
          corresponding DataFrame (or None if input was None/empty). DataFrames
          will have the 'SampleID_Normalized' column added and index reset.
        - A set of unique normalized sample IDs found across all DataFrames.
        - A dictionary mapping normalized sample IDs to the first encountered
          original sample ID.
    """
    dfs = {"genotype": df_genotypes, "alpha": df_pheno_alpha, "numeric": df_pheno_num}
    processed_dfs: Dict[str, Optional[pd.DataFrame]] = {}
    original_id_map: Dict[str, str] = {}  # Map normalized ID -> first original ID
    all_normalized_ids: Set[str] = set()

    for df_type, df in dfs.items():
        if df is not None and not df.empty:
            df_processed = df.copy()
            original_ids = df_processed.index.astype(
                str
            ).tolist()  # Get original IDs before potentially changing index
            df_processed["SampleID_Normalized"] = df_processed.index.map(
                _normalize_sample_id
            )

            all_normalized_ids.update(df_processed["SampleID_Normalized"].unique())

            # Store mapping: normalized -> original
            for i, norm_id in enumerate(df_processed["SampleID_Normalized"]):
                if norm_id not in original_id_map:
                    original_id_map[norm_id] = original_ids[
                        i
                    ]  # Map normalized to the first original ID encountered

            processed_dfs[df_type] = df_processed
        else:
            processed_dfs[df_type] = None

    return processed_dfs, all_normalized_ids, original_id_map


# --- PDF Styling and Content Generation ---


# def _setup_styles() -> Dict[str, ParagraphStyle]:
#     """Creates and returns a dictionary of ReportLab ParagraphStyle objects.

#     Returns:
#         A dictionary mapping style names (e.g., 'body', 'title', 'warning')
#         to their corresponding ParagraphStyle objects.
#     """
#     styles = getSampleStyleSheet()
#     # Reduce body text size slightly for potentially better table fit
#     body_style = ParagraphStyle(
#         name="BodyTextSmall", parent=styles["BodyText"], fontSize=9
#     )
#     table_cell_style = ParagraphStyle(
#         name="TableCell",
#         parent=body_style,
#         leading=11,
#         fontSize=8,  # Smaller font, tighter leading for cells
#     )
#     footer_style = ParagraphStyle(
#         name="FooterStyle",
#         parent=styles["Normal"],
#         alignment=1,
#         fontSize=7,
#         textColor=colors.grey,  # Centered footer
#     )

#     custom_styles = {
#         "body": body_style,
#         "heading": styles["h2"],
#         "title": styles["h1"],
#         "report_title": ParagraphStyle(
#             name="ReportTitle",
#             parent=styles["h1"],
#             alignment=1,
#             fontSize=16,  # Centered H1
#         ),
#         "warning": ParagraphStyle(
#             name="Warning",
#             parent=styles["Heading1"],
#             textColor=colors.red,
#             alignment=1,
#             fontSize=14,
#         ),
#         "table_header": ParagraphStyle(
#             name="TableHeader",
#             parent=styles["BodyText"],
#             fontName="Helvetica-Bold",
#             fontSize=8,  # Smaller header
#         ),
#         "table_cell": table_cell_style,
#         "footer": footer_style,
#     }
#     return custom_styles

def _setup_styles() -> Dict[str, ParagraphStyle]:
    """Creates and returns a dictionary of ReportLab ParagraphStyle objects.

    Returns:
        A dictionary mapping style names (e.g., 'body', 'title', 'warning')
        to their corresponding ParagraphStyle objects.
    """
    styles = getSampleStyleSheet()
    body_style = ParagraphStyle(name='BodyTextSmall', parent=styles['BodyText'], fontSize=9)

    table_cell_style = ParagraphStyle(
        name="TableCell",
        parent=body_style,
        fontSize=7,      
        leading=9,      
    )
    # ****************************************************************

    footer_style = ParagraphStyle(
        name='FooterStyle', parent=styles['Normal'], alignment=1, fontSize=7, textColor=colors.grey # Centered footer
    )

    custom_styles = {
        "body": body_style,
        "heading": styles["h2"],
        "title": styles["h1"],
        "report_title": ParagraphStyle(
            name="ReportTitle", parent=styles["h1"], alignment=1, fontSize=16 # Centered H1
        ),
        "warning": ParagraphStyle(
            name="Warning", parent=styles["Heading1"], textColor=colors.red, alignment=1, fontSize=14
        ),
        "table_header": ParagraphStyle(
            # Keep header size reasonable or reduce slightly if needed too
            name="TableHeader", parent=styles["BodyText"], fontName="Helvetica-Bold", fontSize=8
        ),
        "table_cell": table_cell_style, # Use the modified style
        "footer": footer_style,
    }
    return custom_styles

def _add_page_footer(
    canvas: canvas.Canvas,
    doc: BaseDocTemplate,
    styles: Dict[str, ParagraphStyle],
    UUID: str,
) -> None:
    """Draws the footer on each page."""
    canvas.saveState()
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    footer_text = (
        f"Generated: {timestamp} | Code: v{VERSION} | DB: {DB_VERSION} | UUID: {UUID}"
    )
    p = Paragraph(footer_text, styles["footer"])
    w, h = p.wrap(doc.width, doc.bottomMargin)  # Use doc.width available space
    p.drawOn(
        canvas, doc.leftMargin, doc.bottomMargin - h - 5
    )  # Position below bottom margin
    canvas.restoreState()


def _create_report_header(
    story: List[Flowable], original_id_display: str, styles: Dict[str, ParagraphStyle]
) -> None:
    """Adds standard header elements to the PDF story.

    Args:
        story: The list of ReportLab Flowables to append to.
        original_id_display: The sample ID to display in the header.
        styles: Dictionary containing the required ParagraphStyles.
    """
    story.append(Paragraph("RBCeq2", styles["report_title"]))
    story.append(Spacer(1, 0.1 * inch))
    story.append(Paragraph("Not for clinical use", styles["warning"]))
    story.append(Spacer(1, 0.2 * inch))
    story.append(Paragraph(f"Sample ID: {original_id_display}", styles["heading"]))
    story.append(Spacer(1, 0.1 * inch))  # Reduced space before table title


def _get_data_for_sample(
    norm_id: str, processed_dfs: Dict[str, Optional[pd.DataFrame]]
) -> Tuple[Dict[str, Optional[pd.Series]], Set[str]]:
    """Retrieves data rows for a specific normalized sample ID from processed dfs.

    Args:
        norm_id: The normalized sample ID to retrieve data for.
        processed_dfs: Dictionary mapping data type to processed DataFrame
                       (must contain 'SampleID_Normalized' column).

    Returns:
        A tuple containing:
        - A dictionary mapping data type ('genotype', 'alpha', 'numeric')
          to the pd.Series containing the sample's data (or None if not found).
          The Series index will be the gene/system names.
        - A set of all unique gene/system keys found for this sample across
          all data types.
    """
    sample_data: Dict[str, Optional[pd.Series]] = {}
    all_keys_for_sample: Set[str] = set()

    for df_type, df in processed_dfs.items():
        row_data: Optional[pd.Series] = None
        if df is not None:
            # Find rows matching the normalized ID
            matching_rows = df[df["SampleID_Normalized"] == norm_id]
            if not matching_rows.empty:
                # Select the first match (assuming one row per normalized ID per df type)
                # Drop the helper column to get only gene/system data
                row_data = matching_rows.iloc[0].drop(
                    "SampleID_Normalized", errors="ignore"
                )  # Ignore error if column somehow missing
                all_keys_for_sample.update(row_data.index)

        sample_data[df_type] = row_data

    return sample_data, all_keys_for_sample


def _create_consolidated_table(
    sample_data: Dict[str, Optional[pd.Series]],
    all_keys_for_sample: Set[str],
    styles: Dict[str, ParagraphStyle],
) -> Flowable:
    """Creates the main data table combining genotype and phenotypes.

    Args:
        sample_data: Dict mapping data type to the sample's data Series.
        all_keys_for_sample: Set of all gene/system keys for this sample.
        styles: Dictionary containing the required ParagraphStyles.

    Returns:
        A ReportLab Table object or a Paragraph indicating no data.
    """
    if not all_keys_for_sample:
        return Paragraph("No data available for this sample.", styles["body"])

    # Prepare table header row
    header_content = ["System/Gene", "Genotype", "Phenotype (Alpha)", "Phenotype (Num)"]
    # Wrap header text in Paragraphs using the header style
    header_row: List[Paragraph] = [
        Paragraph(text, styles["table_header"]) for text in header_content
    ]
    combined_table_data: List[List[Flowable]] = [header_row]

    # Sort keys alphabetically for consistent row order
    sorted_keys = sorted(list(all_keys_for_sample))

    for key in sorted_keys:
        # Safely get data using .get(key, default_value) for robustness
        geno_val = (
            sample_data["genotype"].get(key, None)
            if sample_data["genotype"] is not None
            else None
        )
        alpha_val = (
            sample_data["alpha"].get(key, None)
            if sample_data["alpha"] is not None
            else None
        )
        num_val = (
            sample_data["numeric"].get(key, None)
            if sample_data["numeric"] is not None
            else None
        )

        # Format values for the table cells
        # Check data examples to decide on separators if needed
        # Assuming genotype might contain '\n' which _format_cell_content should handle
        formatted_geno = _format_cell_content(geno_val, separator="\n")
        # Use newline as separator if present in data
        formatted_alpha = _format_cell_content(alpha_val)
        # No separator assumed for alpha
        formatted_num = _format_cell_content(
            num_val
        )  # No separator assumed for numeric

        # Create Paragraphs for body cells using the cell style
        combined_table_data.append(
            [
                Paragraph(str(key), styles["table_cell"]),
                Paragraph(formatted_geno, styles["table_cell"]),
                Paragraph(formatted_alpha, styles["table_cell"]),
                Paragraph(formatted_num, styles["table_cell"]),
            ]
        )

    # Slightly adjusted widths
    col_widths = [1.5 * inch, 2.5 * inch, 2.0 * inch, 1.5 * inch]

    # Create the table
    combined_table = Table(combined_table_data, colWidths=col_widths, splitByRow=1)

    # Apply styling (using your original style block)
    combined_table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.grey),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.whitesmoke),
                ("ALIGN", (0, 0), (-1, 0), "CENTER"),  # Header alignment
                ("ALIGN", (0, 1), (0, -1), "LEFT"),  # Gene name left aligned
                ("ALIGN", (1, 1), (-1, -1), "LEFT"),  # Other columns LEFT aligned
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),  # Vertical alignment
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("BOTTOMPADDING", (0, 0), (-1, 0), 8),  # Padding for header
                ("TOPPADDING", (0, 0), (-1, 0), 4),
                ("BOTTOMPADDING", (0, 1), (-1, -1), 4),  # Padding for cells
                ("TOPPADDING", (0, 1), (-1, -1), 4),
                ("BACKGROUND", (0, 1), (-1, -1), colors.beige),
                ("GRID", (0, 0), (-1, -1), 1, colors.black),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
            ]
        )
    )
    return combined_table


def _generate_pdf_report_for_sample(
    norm_id: str,
    original_id_display: str,
    processed_dfs: Dict[str, Optional[pd.DataFrame]],
    styles: Dict[str, ParagraphStyle],
    output_dir: Path,
    UUID: str,
) -> None:
    """Generates and saves a single PDF report for one sample.

    Args:
        norm_id: The normalized sample ID.
        original_id_display: The original sample ID to use for display and filename.
        processed_dfs: Dict mapping data type to the processed DataFrame.
        styles: Dictionary of required ParagraphStyles.
        output_dir: The directory to save the PDF file.
    """
    # Sanitize original ID for filename
    safe_original_id = "".join(
        c if c.isalnum() or c in ("_", "-") else "_" for c in original_id_display
    )
    pdf_filename = os.path.join(output_dir, f"{safe_original_id}_BloodGroupReport.pdf")

    # Use BaseDocTemplate to allow page footers
    doc = BaseDocTemplate(
        pdf_filename,
        pagesize=letter,  # Use letter size as per your original code
        leftMargin=0.75 * inch,
        rightMargin=0.75 * inch,
        topMargin=0.75 * inch,  # Increased top margin slightly
        bottomMargin=0.75 * inch,
    )

    # Define the frame for content flow
    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id="normal")
    # Create a PageTemplate that uses the footer function
    template = PageTemplate(
        id="main",
        frames=[frame],
        onPage=lambda canvas, doc: _add_page_footer(canvas, doc, styles, UUID),
    )
    doc.addPageTemplates([template])

    story: List[Flowable] = []

    # 1. Add Header (after PageTemplate setup)
    _create_report_header(story, original_id_display, styles)

    # 2. Get Data for this sample
    sample_data, all_keys_for_sample = _get_data_for_sample(norm_id, processed_dfs)

    # 3. Add Main Table Title
    story.append(
        Paragraph("Blood Group Genotype & Predicted Phenotype", styles["heading"])
    )
    story.append(Spacer(1, 0.1 * inch))

    # 4. Create and add the data table (or "No data" message)
    table_flowable = _create_consolidated_table(
        sample_data, all_keys_for_sample, styles
    )
    story.append(table_flowable)
    # Removed spacer after table to maximize space

    # 5. Build the PDF
    try:
        doc.build(story)
    except Exception as e:
        logger.error(f"ERROR generating PDF for {norm_id} ({original_id_display}): {e}")


# --- Main Orchestration Function ---


def generate_all_reports(
    df_genotypes: Optional[pd.DataFrame],
    df_pheno_alpha: Optional[pd.DataFrame],
    df_pheno_num: Optional[pd.DataFrame],
    output_name: Path,
    UUID: str,
) -> None:
    """
    Generates consolidated PDF blood group reports for all samples found in input DataFrames.

    Args:
        df_genotypes: DataFrame with genotype data. Index must be SampleID.
        df_pheno_alpha: DataFrame with alphanumeric phenotype data. Index must be SampleID.
        df_pheno_num: DataFrame with numeric phenotype data. Index must be SampleID.
        output_dir: Directory where the generated PDF files will be saved.
    """

    # 1. Prepare Data & Get Unique IDs
    processed_dfs, all_normalized_ids, original_id_map = _prepare_dataframes(
        df_genotypes, df_pheno_alpha, df_pheno_num
    )

    if not all_normalized_ids:
        logger.warning("No samples found in the input DataFrames. Exiting.")
        return

    num_samples = len(all_normalized_ids)
    logger.info(f"Found {num_samples} unique normalized sample IDs.")
    output_dir = f"{output_name}_PDFs"
    # 2. Create output directory
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Reports will be saved to '{os.path.abspath(output_dir)}'.")

    # 3. Setup Styles
    styles = _setup_styles()

    # 4. Iterate Through Each Sample and Generate PDF
    for norm_id in sorted(list(all_normalized_ids)):
        original_id_display = original_id_map.get(
            norm_id, norm_id
        )  # Fallback to norm_id if mapping failed

        try:
            _generate_pdf_report_for_sample(
                norm_id, original_id_display, processed_dfs, styles, output_dir, UUID
            )
        except Exception as e:
            # Catch potential errors during individual PDF generation
            logger.error(
                f"Unhandled exception during PDF generation for {original_id_display}: {e}"
            )
