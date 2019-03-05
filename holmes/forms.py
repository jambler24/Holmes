from django import forms


class UploadFileForm(forms.Form):
    title = forms.CharField(max_length=50)
    file = forms.FileField()


class GeneQueryForm(forms.Form):
    gene_name = forms.CharField(max_length=50)



