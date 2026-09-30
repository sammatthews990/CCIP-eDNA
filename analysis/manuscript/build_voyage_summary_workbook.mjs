import fs from "node:fs/promises";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const outputDir = "outputs/voyage_summary";
const previewDir = "analysis/manuscript/previews";
await fs.mkdir(outputDir, { recursive: true });
await fs.mkdir(previewDir, { recursive: true });

function parseCsv(text) {
  const rows = [];
  let row = [];
  let field = "";
  let quoted = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (quoted) {
      if (ch === '"' && text[i + 1] === '"') {
        field += '"';
        i++;
      } else if (ch === '"') {
        quoted = false;
      } else {
        field += ch;
      }
    } else if (ch === '"') {
      quoted = true;
    } else if (ch === ",") {
      row.push(field);
      field = "";
    } else if (ch === "\n") {
      row.push(field.replace(/\r$/, ""));
      rows.push(row);
      row = [];
      field = "";
    } else {
      field += ch;
    }
  }
  if (field.length || row.length) {
    row.push(field.replace(/\r$/, ""));
    rows.push(row);
  }
  return rows;
}

function excelDate(value) {
  return new Date(`${value}T00:00:00Z`);
}

const summaryRows = parseCsv(await fs.readFile("analysis/manuscript/voyage_summary_table.csv", "utf8"));
const auditRows = parseCsv(await fs.readFile("analysis/manuscript/voyage_summary_audit.csv", "utf8"));

const summaryHeader = [
  "Vessel", "Voyage name/number", "Start date", "End date", "No. reefs visited",
  "No. sites", "Funding program", "Sampling method", "Total samples", "INLA site-model matches"
];
const summaryData = summaryRows.slice(1).map((r) => [
  r[0], r[1], excelDate(r[2]), excelDate(r[3]), Number(r[4]), Number(r[5]),
  r[6], r[7], Number(r[8]), Number(r[9])
]);
const totalSamples = summaryData.reduce((sum, row) => sum + row[8], 0);
const totalMatches = summaryData.reduce((sum, row) => sum + row[9], 0);

const auditHeader = [
  "Voyage", "Observed reef-level designs (design; no. reefs)", "First sample date",
  "Last sample date", "Underlying cull dives", "INLA site-model matches"
];
const auditData = auditRows.slice(1).map((r) => [
  r[0], r[1], excelDate(r[2]), excelDate(r[3]), Number(r[4]), Number(r[5])
]);

const wb = Workbook.create();
const summary = wb.worksheets.add("Voyage summary");
const audit = wb.worksheets.add("Audit");
summary.showGridLines = false;
audit.showGridLines = false;
summary.tabColor = "#0B5563";
audit.tabColor = "#829399";

const font = "Arial";
const dark = "#0B5563";
const accent = "#D97706";
const pale = "#E8F1F3";
const muted = "#5F6B70";

summary.getRange("A2").values = [["Voyages included in the eDNA-COTS analysis"]];
summary.getRange("A2:J2").format.font = { name: font, size: 15, bold: true, color: "#1F2933" };
summary.getRange("A3").values = [[`36 voyages | ${totalSamples.toLocaleString("en-AU")} eDNA samples | ${totalMatches.toLocaleString("en-AU")} INLA site-model matches`]];
summary.getRange("A3:J3").merge();
summary.getRange("A3:J3").format.font = { name: font, size: 10, italic: true, color: muted };
summary.getRange("A4:J4").format.fill = dark;
summary.getRange("A4:J4").format.rowHeight = 3;
summary.getRange("A5:J5").values = [summaryHeader];
summary.getRange("A6:J41").values = summaryData;

const table = summary.tables.add("A5:J41", true, "VoyageSummaryTable");
table.style = "TableStyleMedium2";
table.showBandedRows = true;
table.showFilterButton = true;

summary.getRange("A5:J5").format = {
  fill: dark,
  font: { name: font, size: 10, bold: true, color: "#FFFFFF" },
  horizontalAlignment: "center",
  verticalAlignment: "center",
  wrapText: true
};
summary.getRange("A6:J41").format.font = { name: font, size: 10, color: "#243238" };
summary.getRange("A6:J41").format.verticalAlignment = "center";
summary.getRange("C6:D41").format.numberFormat = "dd mmm yyyy";
summary.getRange("E6:F41").format.numberFormat = "#,##0";
summary.getRange("I6:J41").format.numberFormat = "#,##0";
summary.getRange("C6:F41").format.horizontalAlignment = "right";
summary.getRange("I6:J41").format.horizontalAlignment = "right";

summary.getRange("A42:J42").format = {
  fill: pale,
  font: { name: font, size: 10, bold: true, color: "#1F2933" },
  verticalAlignment: "center"
};
summary.getRange("A42").values = [["Voyage totals"]];
summary.getRange("E42:F42").formulas = [["=SUM(E6:E41)", "=SUM(F6:F41)"]];
summary.getRange("I42:J42").formulas = [["=SUM(I6:I41)", "=SUM(J6:J41)"]];
summary.getRange("E42:F42").format.numberFormat = "#,##0";
summary.getRange("I42:J42").format.numberFormat = "#,##0";
summary.getRange("A44").values = [["Note: reef and site totals are summed voyage-level visits, not unique GBR-wide locations."]];
summary.getRange("A44:J44").merge();
summary.getRange("A44:J44").format.font = { name: font, size: 9, italic: true, color: muted };
summary.getRange("A45").values = [["Source: data/eDNA data_ALL_20260528.xlsx (eDNA_data_ALL and eDNAVoyages sheets)."]];
summary.getRange("A45:J45").merge();
summary.getRange("A45:J45").format.font = { name: font, size: 9, italic: true, color: muted };

summary.getRange("A:A").format.columnWidth = 20;
summary.getRange("B:B").format.columnWidth = 19;
summary.getRange("C:D").format.columnWidth = 13;
summary.getRange("E:F").format.columnWidth = 12;
summary.getRange("G:G").format.columnWidth = 37;
summary.getRange("H:H").format.columnWidth = 34;
summary.getRange("I:I").format.columnWidth = 13;
summary.getRange("J:J").format.columnWidth = 17;
summary.getRange("G6:H41").format.wrapText = true;
summary.getRange("A5:J5").format.rowHeight = 32;
summary.getRange("A6:J41").format.rowHeight = 29;
summary.freezePanes.freezeRows(5);
summary.freezePanes.freezeColumns(2);

for (let i = 0; i < summaryData.length; i++) {
  if (summaryData[i][9] === 0) {
    summary.getRange(`J${i + 6}`).format = {
      fill: "#F3F4F6",
      font: { name: font, size: 10, color: muted },
      horizontalAlignment: "right"
    };
  }
}

audit.getRange("A2").values = [["Voyage table audit and definitions"]];
audit.getRange("A2:F2").format.font = { name: font, size: 15, bold: true, color: "#1F2933" };
audit.getRange("A3:F3").format.fill = dark;
audit.getRange("A3:F3").format.rowHeight = 3;
audit.getRange("A4").values = [["Sampling design is classified from the observed number of sites and sample rows per site within each reef."]];
audit.getRange("A5").values = [["INLA site-model matches use one unique eDNA campaign per first subsequent cull-site visit, within 2,000 m and 183 days."]];
audit.getRange("A6").values = [["Each of the 519 site visits appears once. The audit also reports the 1,064 underlying cull dives aggregated into those responses."]];
audit.getRange("A4:F6").format.font = { name: font, size: 10, italic: true, color: muted };
audit.getRange("A8:F8").values = [auditHeader];
audit.getRange("A9:F44").values = auditData;
const auditTable = audit.tables.add("A8:F44", true, "VoyageAuditTable");
auditTable.style = "TableStyleMedium2";
auditTable.showBandedRows = true;
auditTable.showFilterButton = true;
audit.getRange("A8:F8").format = {
  fill: dark,
  font: { name: font, size: 10, bold: true, color: "#FFFFFF" },
  horizontalAlignment: "center",
  verticalAlignment: "center",
  wrapText: true
};
audit.getRange("A9:F44").format.font = { name: font, size: 10, color: "#243238" };
audit.getRange("A9:F44").format.verticalAlignment = "center";
audit.getRange("C9:D44").format.numberFormat = "dd mmm yyyy";
audit.getRange("E9:F44").format.numberFormat = "#,##0";
audit.getRange("C9:F44").format.horizontalAlignment = "right";
audit.getRange("A46").values = [["Sources"]];
audit.getRange("A46").format.font = { name: font, size: 10, bold: true, color: dark };
audit.getRange("A47").values = [["eDNA data and voyage metadata: data/eDNA data_ALL_20260528.xlsx"]];
audit.getRange("A48").values = [["Cull records: data/260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx, Cull sheet"]];
audit.getRange("A47:F48").format.font = { name: font, size: 9, italic: true, color: muted };
audit.getRange("A:A").format.columnWidth = 20;
audit.getRange("B:B").format.columnWidth = 58;
audit.getRange("C:D").format.columnWidth = 16;
audit.getRange("E:F").format.columnWidth = 17;
audit.getRange("A8:F8").format.rowHeight = 32;
audit.getRange("A9:F44").format.rowHeight = 24;
audit.freezePanes.freezeRows(8);
audit.freezePanes.freezeColumns(1);

wb.recalculate();

const summaryInspect = await wb.inspect({
  kind: "table", range: "Voyage summary!A2:J45", include: "values,formulas",
  tableMaxRows: 45, tableMaxCols: 10
});
console.log(summaryInspect.ndjson);
const errorInspect = await wb.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!",
  options: { useRegex: true, maxResults: 300 },
  summary: "final formula error scan"
});
console.log(errorInspect.ndjson);

for (const sheetName of ["Voyage summary", "Audit"]) {
  const preview = await wb.render({ sheetName, autoCrop: "all", scale: 1, format: "png" });
  const safeName = sheetName.toLowerCase().replaceAll(" ", "_");
  await fs.writeFile(`${previewDir}/${safeName}.png`, new Uint8Array(await preview.arrayBuffer()));
}

const out = await SpreadsheetFile.exportXlsx(wb);
await out.save(`${outputDir}/voyage_summary_for_manuscript_site_model_aligned.xlsx`);
console.log(`${outputDir}/voyage_summary_for_manuscript_site_model_aligned.xlsx`);
