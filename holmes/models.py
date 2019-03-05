from __future__ import unicode_literals

from django.db import models

# Create your models here.


class Experiment(models.Model):
	"""This describes the experimental setup of the gene expression data"""

	# Fields

	experiment_shortname = models.CharField(max_length=30, help_text='Enter field documentation')
	comparisons = models.TextField(help_text='Enter the conditions separated by a comma.')
	time_series = models.BooleanField()

	def __str__(self):
		return self.experiment_shortname


class CurrentSettings(models.Model):
	"""The settings for the current user and session"""

	# Fields

	current_experimental_setup = models.CharField(max_length=30, help_text='Currently selected experimental conditions')
