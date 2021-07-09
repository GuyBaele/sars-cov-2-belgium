import matplotlib.pyplot as plt  # type: ignore
import pandas as pd  # type: ignore

df = pd.read_csv("data/metadata.tsv", sep="\t")

df = df[df["country"] == "India"]
df = df[df["pango_lineage"] == "B.1.617.2"]
df = df[df["date"].apply(lambda x: len(x) == 10)]
df = df[df["date"].apply(lambda x: "X" not in str(x))]
df["date"] = df["date"].astype("datetime64")

for index, row in df[df["date"].dt.isocalendar().week == 53].iterrows():
    print(row)

plt.figure(figsize=(20, 10))
ax = (
    df["strain"]
    .groupby([df["date"].dt.isocalendar().year, df["date"].dt.isocalendar().week])
    .count()
    .plot(kind="bar")
)
ax.set_facecolor("#eeeeee")
ax.set_xlabel("month")
ax.set_ylabel("count")
ax.set_title("Indian Delta genomes submitted to GISAID")
plt.show()
