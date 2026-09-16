import hashlib
import html
import os
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

# --- Reference Alleles ---
# Used to determine if a genotype is homozygous reference (ref/ref).
# If both alleles in a genotype match the reference for that system,
# the row is NOT highlighted. Any deviation is subtly highlighted.

REFERENCE_ALLELES: Dict[str, str] = {
    "CROM": "CROM*01",
    "FY": "FY*01",
    "KN": "KN*01",
    "RHD": "RHD*01",
    "SC": "SC*01",
    "VEL": "VEL*01",
    "CD59": "CD59*01",
    "IN": "IN*02",
    "RAPH": "RAPH*01",
    "DO": "DO*02",
    "ABCC4": "ABCC4*01",
    "JMH": "JMH*01",
    "ER": "ER*01",
    "DI": "DI*02",
    "SID": "SID*01",
    "JK": "JK*01",
    "CTL2": "CTL2*01",
    "FUT1": "FUT1*01",
    "FUT2": "FUT2*01",
    "FUT3": "FUT3*01.01",
    "KLF1": "KLF1*01",
    "LU": "LU*02",
    "LW": "LW*05",
    "MAM": "MAM*01",
    "OK": "OK*01.01",
    "ABCB6": "ABCB6*01",
    "GE": "GE*01",
    "KANNO": "KANNO*01",
    "A4GALT": "A4GALT*01",
    "GLOB": "GLOB*01",
    "ABCG2": "ABCG2*01",
    "PIGG": "PIGG*01",
    "AUG": "AUG*01",
    "C4A": "C4A*03",
    "C4B": "C4B*03",
    "GCNT2": "GCNT2*01",
    "RHAG": "RHAG*01",
    "CO": "CO*01.01",
    "KEL": "KEL*02",
    "YT": "YT*01",
    "GBGT1": "GBGT1*01N.01",
    "GIL": "GIL*01",
    "GATA1": "GATA1*01",
    "XG": "XG*01",
    "CD99": "CD99*01",
    "XK": "XK*01",
    "ABO": "ABO*A1.01",
    "ABCC1": "ABCC1*01",
    "GYPA": "GYPA*01",
    "GYPB": "GYPB*04",
    "GYP": "GYP*00",
    "RHCE": "RHCE*01",
    "ATP11C": "ATP11C*01",
    # HPA reference alleles
    "HPA1": "HPA1*01",
    "HPA2": "HPA2*01",
    "HPA3": "HPA3*01",
    "HPA4": "HPA4*01",
    "HPA5": "HPA5*01",
    "HPA6": "HPA6*01",
    "HPA7": "HPA7*01",
    "HPA8": "HPA8*01",
    "HPA9": "HPA9*01",
    "HPA10": "HPA10*01",
    "HPA11": "HPA11*01",
    "HPA12": "HPA12*01",
    "HPA13": "HPA13*01",
    "HPA14": "HPA14*01",
    "HPA15": "HPA15*01",
    "HPA16": "HPA16*01",
    "HPA17": "HPA17*01",
    "HPA18": "HPA18*01",
    "HPA19": "HPA19*01",
    "HPA20": "HPA20*01",
    "HPA21": "HPA21*01",
    "HPA22": "HPA22*01",
    "HPA23": "HPA23*01",
    "HPA24": "HPA24*01",
    "HPA25": "HPA25*01",
    "HPA26": "HPA26*01",
    "HPA27": "HPA27*01",
    "HPA28": "HPA28*01",
    "HPA29": "HPA29*01",
    "HPA30": "HPA30*01",
    "HPA31": "HPA31*01",
    "HPA32": "HPA32*01",
    "HPA33": "HPA33*01",
    "HPA34": "HPA34*01",
    "HPA35": "HPA35*01",
}

# Highlight color for non-reference genotypes — visible yellow
HIGHLIGHT_COLOR = colors.Color(1.0, 0.92, 0.40)  # Bright yellow
# Default row backgrounds for alternating stripes
DEFAULT_ROW_COLOR_EVEN = colors.Color(0.95, 0.95, 0.93)  # Beige-ish
DEFAULT_ROW_COLOR_ODD = colors.Color(1.0, 1.0, 1.0)  # White


# --- Helper Functions ---


def _normalize_sample_id(sample_id: str) -> str:
    """Validate an exact sample ID without removing suffixes or converting types.

    Raises:
        ValueError: The identity is empty or is not a string.
    """
    if not isinstance(sample_id, str) or not sample_id:
        raise ValueError(f"PDF sample IDs must be nonempty strings: {sample_id!r}")
    return sample_id


def _pdf_filename(sample_id: str) -> str:
    """Build a bounded ASCII filename tied to the complete, exact sample ID.

    The readable prefix is for navigation; the full digest distinguishes IDs
    whose prefixes sanitize alike, differ only in case, or exceed its length.
    The name depends only on the identity, not the other samples in a run.
    """
    sample_id = _normalize_sample_id(sample_id)
    prefix = "".join(
        c if c.isascii() and (c.isalnum() or c in "_-") else "_"
        for c in sample_id[:48]
    )
    digest = hashlib.sha256(sample_id.encode("utf-8")).hexdigest()
    return f"{prefix}__{digest}_BloodGroupReport.pdf"


def _format_cell_content(text: Optional[Any], separator: Optional[str] = None) -> str:
    """Formats cell text for PDF tables, handling N/A, line breaks, and HTML escaping."""
    if pd.isna(text) or text == "" or text is None:
        return "N/A"

    text_str = str(text)
    escaped_text = html.escape(text_str)

    if separator:
        items = text_str.split(separator)
        escaped_items = [html.escape(item.strip()) for item in items]
        return "<br/>".join(escaped_items)
    else:
        return escaped_text


def _is_hpa_key(key: str) -> bool:
    """Returns True if the gene/system key is an HPA allele."""
    return key.upper().startswith("HPA")


def _is_ref_ref(key: str, genotype_value: Optional[Any]) -> bool:
    """Checks if a genotype value represents a homozygous reference call.

    Parses the genotype string (expected format: "ALLELE1/ALLELE2") and checks
    whether both alleles match the known reference allele for this system/gene.

    Args:
        key: The gene/system name (e.g., "FY", "KEL", "HPA1").
        genotype_value: The raw genotype value from the DataFrame.

    Returns:
        True if the genotype is homozygous reference, False otherwise.
    """
    if pd.isna(genotype_value) or genotype_value is None or str(genotype_value).strip() == "":
        return True  # Treat missing/empty as "nothing to highlight"

    geno_str = str(genotype_value).strip()
    ref_allele = REFERENCE_ALLELES.get(key)

    if ref_allele is None:
        # Unknown system — don't highlight if we don't know the reference
        return True

    # Expected format: "ALLELE1/ALLELE2" but may also contain commas
    # (e.g. multi-allele calls like "ALLELE1/ALLELE2,ALLELE3/ALLELE4").
    # Remove any leading/trailing whitespace and newlines for the check.
    clean_geno = geno_str.replace("\n", "").replace("\r", "").strip()

    # If comma-separated (multiple loci), check each pair independently
    segments = [s.strip() for s in clean_geno.split(",")]
    for segment in segments:
        alleles = [a.strip() for a in segment.split("/")]
        if len(alleles) == 2:
            if not (alleles[0] == ref_allele and alleles[1] == ref_allele):
                return False
        elif len(alleles) == 1 and alleles[0]:
            # Single allele listed — check if it matches ref
            if alleles[0] != ref_allele:
                return False
        # Skip empty segments

    return True


# --- Data Preparation ---


def _prepare_dataframes(
    df_genotypes: Optional[pd.DataFrame],
    df_pheno_alpha: Optional[pd.DataFrame],
    df_pheno_num: Optional[pd.DataFrame],
) -> Tuple[Dict[str, Optional[pd.DataFrame]], Set[str], Dict[str, str]]:
    """Copy report frames while preserving exact, unique string sample identities.

    No metadata column is inserted into the biological data. Missing rows in a
    phenotype frame stay missing rather than borrowing another sample's row.

    Raises:
        ValueError: A frame contains an invalid or duplicate sample identity.
    """
    dfs = {"genotype": df_genotypes, "alpha": df_pheno_alpha, "numeric": df_pheno_num}
    processed_dfs: Dict[str, Optional[pd.DataFrame]] = {}
    original_id_map: Dict[str, str] = {}
    all_sample_ids: Set[str] = set()

    for df_type, df in dfs.items():
        if df is not None and not df.empty:
            sample_ids = [_normalize_sample_id(sample) for sample in df.index]
            if not df.index.is_unique:
                duplicates = sorted(set(df.index[df.index.duplicated()]))
                raise ValueError(f"Duplicate PDF sample IDs in {df_type}: {duplicates!r}")
            processed_dfs[df_type] = df.copy()
            all_sample_ids.update(sample_ids)
            original_id_map.update({sample: sample for sample in sample_ids})
        else:
            processed_dfs[df_type] = None

    return processed_dfs, all_sample_ids, original_id_map


# --- PDF Styling and Content Generation ---


def _setup_styles() -> Dict[str, ParagraphStyle]:
    """Creates and returns a dictionary of ReportLab ParagraphStyle objects."""
    styles = getSampleStyleSheet()
    body_style = ParagraphStyle(
        name="BodyTextSmall", parent=styles["BodyText"], fontSize=9
    )
    table_cell_style = ParagraphStyle(
        name="TableCell",
        parent=body_style,
        fontSize=7,
        leading=9,
    )
    footer_style = ParagraphStyle(
        name="FooterStyle",
        parent=styles["Normal"],
        alignment=1,
        fontSize=7,
        textColor=colors.grey,
    )
    section_heading_style = ParagraphStyle(
        name="SectionHeading",
        parent=styles["h3"],
        fontSize=11,
        spaceAfter=4,
        spaceBefore=10,
    )

    summary_style = ParagraphStyle(
        name="SummaryStyle",
        parent=body_style,
        fontSize=9,
        textColor=colors.Color(0.3, 0.3, 0.3),
        spaceBefore=2,
        spaceAfter=6,
    )

    custom_styles = {
        "body": body_style,
        "heading": styles["h2"],
        "title": styles["h1"],
        "report_title": ParagraphStyle(
            name="ReportTitle",
            parent=styles["h1"],
            alignment=1,
            fontSize=16,
        ),
        "warning": ParagraphStyle(
            name="Warning",
            parent=styles["Heading1"],
            textColor=colors.red,
            alignment=1,
            fontSize=14,
        ),
        "table_header": ParagraphStyle(
            name="TableHeader",
            parent=styles["BodyText"],
            fontName="Helvetica-Bold",
            fontSize=8,
        ),
        "table_cell": table_cell_style,
        "footer": footer_style,
        "section_heading": section_heading_style,
        "summary": summary_style,
    }
    return custom_styles


def _add_page_footer(
    canvas: canvas.Canvas,
    doc: BaseDocTemplate,
    styles: Dict[str, ParagraphStyle],
    UUID: str,
) -> None:
    """Draws the footer on each page, including page number."""
    canvas.saveState()
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    page_num = canvas.getPageNumber()
    footer_text = (
        f"Generated: {timestamp} | Code: v{VERSION} | DB: {DB_VERSION} | UUID: {UUID}"
    )
    p = Paragraph(footer_text, styles["footer"])
    w, h = p.wrap(doc.width, doc.bottomMargin)
    p.drawOn(canvas, doc.leftMargin, doc.bottomMargin - h - 5)

    # Page number on the right side
    canvas.setFont("Helvetica", 7)
    canvas.setFillColor(colors.grey)
    canvas.drawRightString(
        doc.width + doc.leftMargin,
        doc.bottomMargin - h - 5,
        f"Page {page_num}",
    )
    canvas.restoreState()


def _create_report_header(
    story: List[Flowable], original_id_display: str, styles: Dict[str, ParagraphStyle]
) -> None:
    """Adds standard header elements to the PDF story."""
    story.append(Paragraph("RBCeq2", styles["report_title"]))
    story.append(Spacer(1, 0.1 * inch))
    story.append(Paragraph("Not for clinical use", styles["warning"]))
    story.append(Spacer(1, 0.2 * inch))
    story.append(Paragraph(f"Sample ID: {html.escape(original_id_display)}", styles["heading"]))
    story.append(Spacer(1, 0.1 * inch))


def _get_data_for_sample(
    norm_id: str, processed_dfs: Dict[str, Optional[pd.DataFrame]]
) -> Tuple[Dict[str, Optional[pd.Series]], Set[str]]:
    """Retrieve rows for an exact sample ID from the validated report frames."""
    sample_data: Dict[str, Optional[pd.Series]] = {}
    all_keys_for_sample: Set[str] = set()

    for df_type, df in processed_dfs.items():
        row_data: Optional[pd.Series] = None
        if df is not None and norm_id in df.index:
            row_data = df.loc[norm_id]
            all_keys_for_sample.update(row_data.index)

        sample_data[df_type] = row_data

    return sample_data, all_keys_for_sample


def _build_table_from_keys(
    keys: List[str],
    sample_data: Dict[str, Optional[pd.Series]],
    styles: Dict[str, ParagraphStyle],
) -> Tuple[Table, int]:
    """Builds a ReportLab Table from a list of gene/system keys.

    Returns the Table and the count of non-ref/ref rows (highlighted rows).
    """
    header_content = ["System/Gene", "Genotype", "Phenotype (Alpha)", "Phenotype (Num)"]
    header_row: List[Paragraph] = [
        Paragraph(text, styles["table_header"]) for text in header_content
    ]
    table_data: List[List[Flowable]] = [header_row]
    highlight_rows: List[int] = []  # 1-based row indices to highlight

    for i, key in enumerate(keys):
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

        formatted_geno = _format_cell_content(geno_val, separator=",")
        formatted_alpha = _format_cell_content(alpha_val, separator="|")
        formatted_num = _format_cell_content(num_val, separator="|")

        # Check if this row should be highlighted (not ref/ref)
        if not _is_ref_ref(key, geno_val):
            highlight_rows.append(i + 1)  # +1 because row 0 is the header

        table_data.append(
            [
                Paragraph(str(key), styles["table_cell"]),
                Paragraph(formatted_geno, styles["table_cell"]),
                Paragraph(formatted_alpha, styles["table_cell"]),
                Paragraph(formatted_num, styles["table_cell"]),
            ]
        )

    col_widths = [1.5 * inch, 2.5 * inch, 2.0 * inch, 1.5 * inch]
    table = Table(table_data, colWidths=col_widths, splitByRow=1)

    # Base table style
    style_commands = [
        ("BACKGROUND", (0, 0), (-1, 0), colors.grey),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.whitesmoke),
        ("ALIGN", (0, 0), (-1, 0), "CENTER"),
        ("ALIGN", (0, 1), (0, -1), "LEFT"),
        ("ALIGN", (1, 1), (-1, -1), "LEFT"),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("BOTTOMPADDING", (0, 0), (-1, 0), 8),
        ("TOPPADDING", (0, 0), (-1, 0), 4),
        ("BOTTOMPADDING", (0, 1), (-1, -1), 4),
        ("TOPPADDING", (0, 1), (-1, -1), 4),
        ("GRID", (0, 0), (-1, -1), 1, colors.black),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
    ]

    # Alternating row backgrounds (applied first, then overridden by highlights)
    for row_idx in range(1, len(table_data)):
        bg = DEFAULT_ROW_COLOR_EVEN if (row_idx % 2 == 1) else DEFAULT_ROW_COLOR_ODD
        style_commands.append(
            ("BACKGROUND", (0, row_idx), (-1, row_idx), bg)
        )

    # Apply yellow highlight to non-ref/ref rows (overrides alternating colour)
    for row_idx in highlight_rows:
        style_commands.append(
            ("BACKGROUND", (0, row_idx), (-1, row_idx), HIGHLIGHT_COLOR)
        )

    table.setStyle(TableStyle(style_commands))
    return table, len(highlight_rows)


def _create_consolidated_table(
    sample_data: Dict[str, Optional[pd.Series]],
    all_keys_for_sample: Set[str],
    styles: Dict[str, ParagraphStyle],
) -> Tuple[List[Flowable], int, int]:
    """Creates the main data tables, separating HPA alleles to the bottom.

    Returns a tuple of:
    - List of Flowables (blood group table, then optionally HPA section)
    - Count of non-ref blood group systems
    - Count of non-ref HPA systems
    """
    if not all_keys_for_sample:
        return [Paragraph("No data available for this sample.", styles["body"])], 0, 0

    sorted_keys = sorted(list(all_keys_for_sample))

    # Separate HPA keys from blood group keys
    blood_group_keys = [k for k in sorted_keys if not _is_hpa_key(k)]
    hpa_keys = [k for k in sorted_keys if _is_hpa_key(k)]

    flowables: List[Flowable] = []
    bg_highlight_count = 0
    hpa_highlight_count = 0

    # Blood group systems table
    if blood_group_keys:
        bg_table, bg_highlight_count = _build_table_from_keys(
            blood_group_keys, sample_data, styles
        )
        flowables.append(bg_table)

    # HPA table (separated at the bottom)
    if hpa_keys:
        flowables.append(Spacer(1, 0.2 * inch))
        flowables.append(
            Paragraph("Human Platelet Antigens (HPA)", styles["section_heading"])
        )
        flowables.append(Spacer(1, 0.05 * inch))
        hpa_table, hpa_highlight_count = _build_table_from_keys(
            hpa_keys, sample_data, styles
        )
        flowables.append(hpa_table)

    return flowables, bg_highlight_count, hpa_highlight_count


def _generate_pdf_report_for_sample(
    norm_id: str,
    original_id_display: str,
    processed_dfs: Dict[str, Optional[pd.DataFrame]],
    styles: Dict[str, ParagraphStyle],
    output_dir: Path,
    UUID: str,
) -> None:
    """Generate one PDF report, propagating rendering and write failures."""
    pdf_filename = os.path.join(output_dir, _pdf_filename(original_id_display))

    doc = BaseDocTemplate(
        pdf_filename,
        pagesize=letter,
        leftMargin=0.75 * inch,
        rightMargin=0.75 * inch,
        topMargin=0.75 * inch,
        bottomMargin=0.75 * inch,
    )

    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id="normal")
    template = PageTemplate(
        id="main",
        frames=[frame],
        onPage=lambda canvas, doc: _add_page_footer(canvas, doc, styles, UUID),
    )
    doc.addPageTemplates([template])

    story: List[Flowable] = []

    # 1. Header
    _create_report_header(story, original_id_display, styles)

    # 2. Get data for this sample
    sample_data, all_keys_for_sample = _get_data_for_sample(norm_id, processed_dfs)

    # 3. Main table title
    story.append(
        Paragraph("Blood Group Genotype & Predicted Phenotype", styles["heading"])
    )
    story.append(Spacer(1, 0.1 * inch))

    # 4. Colour key / legend
    legend_data = [
        [
            Paragraph("", styles["table_cell"]),
            Paragraph("Homozygous reference genotype", styles["table_cell"]),
        ],
        [
            Paragraph("", styles["table_cell"]),
            Paragraph(
                "Non-reference genotype detected", styles["table_cell"]
            ),
        ],
    ]
    legend_table = Table(legend_data, colWidths=[0.3 * inch, 3.0 * inch])
    legend_table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (0, 0), DEFAULT_ROW_COLOR_EVEN),
                ("BACKGROUND", (0, 1), (0, 1), HIGHLIGHT_COLOR),
                ("BOX", (0, 0), (0, 0), 0.75, colors.black),
                ("BOX", (0, 1), (0, 1), 0.75, colors.black),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                ("TOPPADDING", (0, 0), (-1, -1), 2),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
            ]
        )
    )
    story.append(legend_table)
    story.append(Spacer(1, 0.15 * inch))

    # 5. Create tables (blood groups + HPA separated)
    table_flowables, bg_count, hpa_count = _create_consolidated_table(
        sample_data, all_keys_for_sample, styles
    )

    # 5a. Summary count of non-reference genotypes
    sorted_keys = sorted(list(all_keys_for_sample))
    total_bg = len([k for k in sorted_keys if not _is_hpa_key(k)])
    total_hpa = len([k for k in sorted_keys if _is_hpa_key(k)])

    summary_parts = []
    if total_bg > 0:
        summary_parts.append(
            f"{bg_count} of {total_bg} blood group system{'s' if total_bg != 1 else ''}"
            f" show{'s' if bg_count == 1 else ''} non-reference genotype{'s' if bg_count != 1 else ''}"
        )
    if total_hpa > 0:
        summary_parts.append(
            f"{hpa_count} of {total_hpa} HPA system{'s' if total_hpa != 1 else ''}"
            f" show{'s' if hpa_count == 1 else ''} non-reference genotype{'s' if hpa_count != 1 else ''}"
        )
    if summary_parts:
        summary_text = " | ".join(summary_parts)
        story.append(Paragraph(summary_text, styles["summary"]))
        story.append(Spacer(1, 0.05 * inch))

    story.extend(table_flowables)

    # 6. Build the PDF
    doc.build(story)


# --- Main Orchestration Function ---


def generate_all_reports(
    df_genotypes: Optional[pd.DataFrame],
    df_pheno_alpha: Optional[pd.DataFrame],
    df_pheno_num: Optional[pd.DataFrame],
    output_name: Path,
    UUID: str,
) -> None:
    """Generate every sample report and report any failures to the caller.

    Successful reports are retained when another sample fails. Each failure is
    logged with its sample identity before a combined error is raised.

    Raises:
        ValueError: Sample identities are invalid, duplicated, or map to one filename.
        RuntimeError: One or more sample reports could not be generated.
    """
    processed_dfs, all_sample_ids, original_id_map = _prepare_dataframes(
        df_genotypes, df_pheno_alpha, df_pheno_num
    )

    if not all_sample_ids:
        logger.warning("No samples found in the input DataFrames. Exiting.")
        return

    num_samples = len(all_sample_ids)
    logger.info(f"Found {num_samples} unique sample IDs.")
    output_dir = f"{output_name}_PDFs"

    filenames: Dict[str, str] = {}
    for sample_id in sorted(all_sample_ids):
        name = _pdf_filename(sample_id).casefold()
        if name in filenames:
            raise ValueError(
                f"PDF filename collision for {filenames[name]!r} and {sample_id!r}"
            )
        filenames[name] = sample_id

    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Reports will be saved to '{os.path.abspath(output_dir)}'.")

    styles = _setup_styles()

    failed_samples: List[str] = []
    for norm_id in sorted(list(all_sample_ids)):
        original_id_display = original_id_map.get(norm_id, norm_id)

        try:
            _generate_pdf_report_for_sample(
                norm_id, original_id_display, processed_dfs, styles, output_dir, UUID
            )
        except Exception as e:
            logger.error(
                f"PDF generation failed for {original_id_display!r}: {e}"
            )
            failed_samples.append(original_id_display)

    if failed_samples:
        names = ", ".join(repr(sample) for sample in failed_samples)
        raise RuntimeError(
            f"Failed to generate {len(failed_samples)} PDF report(s): {names}. "
            "See the log for each failure."
        )
